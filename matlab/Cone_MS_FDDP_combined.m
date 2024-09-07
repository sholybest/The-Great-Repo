%%
clear; clc;

%%
global N; global dt_mpc; global m_mpc; global IB IB_INV; 
global init_leg_position leg_position I3x12 I3; global g; 
global mu_friction; global fz_max; global fz_min;
global weight_p weight_v weight_psi weight_w weight_u weight_cone;
global terminal_weight_p terminal_weight_v terminal_weight_psi terminal_weight_w;
global backtracking_alpha step_size_alpha;
global p_ref v_ref R_ref w_ref;
global nom_error_p nom_error_v nom_error_psi nom_error_w nominal_p nominal_v nominal_R nominal_w;
global p0 v0 R0 w0;
global state_p state_v state_R state_w;
global regularization_mu;
global state_gap state_gap_init delta1 delta2;
global Qu_data Quu_data Vx_data Vxx_data Qux_data;
global isProjected cone_feasibility; 

%% Parameters
N = 50;
dt_mpc = 0.02;
m_mpc = 40;
IB = diag([0.4 2.1 2.1]);
IB_INV = inv(IB);
mu_friction = 0.6; fz_max = 666; fz_min = 10;
I3 = eye(3);

isProjected = [0;0;0;0];
Qu_data = zeros(12,N);
Quu_data = zeros(12,12*N);
Vx_data = zeros(12,N);
Vxx_data = zeros(12,12*N);
state_gap = zeros(12,N+1);
delta1 = 0; delta2 = 0;
cone_feasibility = 1;
step_size_alpha = 0.8;
regularization_mu = 1e-9;
MAX_ITER = 500;
MAX_SUB_ITER = 10;
epsilon_main = 1e-3;
epsilon_gap = 1e-9;

weight_p   = [1;1;50]; % 1;1;50?
weight_v   = [1;1;1];  
weight_psi = [1;1;1]; % 5;5;5 or 60;60;60
weight_w   = [1e-2;1e-2;1e-2];
weight_u   = [1;1;1; 1;1;1; 1;1;1; 1;1;1]*1e-6;

terminal_weight_p   = weight_p;
terminal_weight_v   = weight_v;
terminal_weight_psi = weight_psi;
terminal_weight_w   = weight_w;

%% Initialization
r1 = [ 0.317;  0.12; -0.48];
r2 = [ 0.317; -0.12; -0.48];
r3 = [-0.317;  0.12; -0.48];
r4 = [-0.317; -0.12; -0.48];

% r1 = [ 0.;  0.12; -0.48];
% r2 = [ 0.; -0.12; -0.48];
% r3 = [-0.;  0.12; -0.48];
% r4 = [-0.; -0.12; -0.48];

init_leg_position = [r1; r2; r3; r4]; % 발 위치가 변할 경우 옆으로 이어붙이기.

I3x12 = repmat(eye(3),1,4);
g = [0;0;-9.81];

p_ref = repmat([0;0;0.5],1,N);
v_ref = repmat([0;0;0],1,N);
R_ref = repmat(eye(3),1,N); 
w_ref = repmat([0;0;0],1,N);

p0 = [0;0;0.48];
v0 = [0.01;-0.1;-0.01];
R0 = rotz(0.1)*roty(-.8)*rotx(.1);
w0 = [0.01;-0.01;0.01];

% p0 = [0;0.1;0.27];
% v0 = [0.05;-0.1;-0.04];
% R0 = rotz(15)*roty(12)*rotx(14);
% w0 = [0.1;0.1;0.1];

state_p = zeros(3,N+1);
state_v = zeros(3,N+1);
state_R = zeros(3,3*(N+1)); % 3x3N
state_w = zeros(3,N+1);

state_gap = zeros(12,N+1);
state_gap_init = zeros(12,N+1);

error_p   = zeros(3,N);
error_v   = zeros(3,N);
error_psi = zeros(3,N);
error_w   = zeros(3,N);

%% DDP iterations
% initial nominal trajectory
initial_input_U = repmat([0;0;100; 0;0;100; 0;0;100; 0;0;100],1,N); 
% initial_input_U = repmat([0;0;100; 0;0;0; 0;0;0; 0;0;100],1,N);
contact_sequence = repmat([1;1;1;1],1,N);
% contact_sequence = repmat([1;0;0;1],1,N);

nominal_p = repmat(p0,1,N);
nominal_v = repmat(v0,1,N);
nominal_R = repmat(R0,1,N);
nominal_w = repmat(w0,1,N);

state_p = [p0, nominal_p];
state_v = [v0, nominal_v];
state_R = [R0, nominal_R];
state_w = [w0, nominal_w];

%%
objective = Objective_Function(state_p, state_v, state_R, state_w, initial_input_U); obj_data(1) = objective;

tic;

prev_objective = 0;
initial_objective = objective;
minimum_objective = objective;
input_U = initial_input_U;
optimal_U = initial_input_U+1;
iter = 0;

while 1
    iter = iter + 1; fprintf("\n\niteration "+string(iter)+":\n");    
    
    % Gap calculations
    for k = 2:N+1
        [rollout_p, rollout_v, rollout_R, rollout_w] = Euler_Integration([state_p(:,k-1),state_v(:,k-1),state_R(:,3*(k-1)-2:3*(k-1)),state_w(:,k-1)], input_U(:,k-1));
        state_gap_init(:,k) = [rollout_p-state_p(:,k);
                          rollout_v-state_v(:,k);
                          Log_mat( state_R(:,3*k-2:3*k).'*rollout_R );
                          rollout_w - (state_R(:,3*k-2:3*k).'*rollout_R).'*state_w(:,k);
                         ];
    end
    state_gap = state_gap_init;
    
    % Feasibility_Check(input_U, contact_sequence);
    cone_feasibility = mod(iter,2) == 0
    if cone_feasibility
        step_size_alpha = 1.0;
    else
        step_size_alpha = 0.5;
    end
    
    expected_improvement = delta1*step_size_alpha + 0.5*delta2*step_size_alpha

    % Backward pass
    [del_u_feedforward, del_u_feedback_gain] = BackPass_DDP(input_U, contact_sequence);

    % Forward Pass
    [optimal_U, optimal_p, optimal_v, optimal_R, optimal_w] = ForwardPass_DDP(del_u_feedforward, del_u_feedback_gain, input_U, contact_sequence);
    objective_forwardpass = Objective_Function([p0,optimal_p], [v0,optimal_v], [R0, optimal_R], [w0,optimal_w], optimal_U);

    % Next step update
    input_U = optimal_U;
    state_p(:,2:end) = optimal_p;
    state_v(:,2:end) = optimal_v;
    state_R(:,4:end) = optimal_R;
    state_w(:,2:end) = optimal_w;

    prev_objective = objective;
    objective = objective_forwardpass; obj_data(iter+1) = objective;

    % exit conditions
    residual = ((prev_objective-objective)/initial_objective)^2;
    fprintf("residual: "+string(residual)+"\n");
    fprintf("max gap: "+string(max(max(abs(state_gap))))+"\n");
    exit_condition1 = (max(max(abs(state_gap))) < epsilon_gap) && residual<epsilon_main;
    exit_condition2 = iter >= MAX_ITER;

    if ( exit_condition1 | exit_condition2 )
        break;
    end

end
fprintf("\n");
toc;

%%
for k=1:N
    error_psi(:,k) = Log_mat( R_ref(:,3*k-2:3*k).' * state_R(:,3*(k+1)-2:3*(k+1)) );
end

%% plots
figure; plot(0:iter,obj_data); grid;
title("Objective Function"); xlabel("iterations"); ylabel("Objective");

figure; plot(vecnorm(error_psi*180/pi)); grid; 
xlim([1 N]); %ylim([0 25]);
title("Rotational error"); xlabel("steps"); ylabel("rotational error (deg)");

figure; hold on; grid; xlim([1,N+1]);
plot(state_p(1,:),"r"); plot(state_p(2,:),"g"); plot(state_p(3,:),"b");
title("Position"); legend("px","py","pz","Location","northwest");

figure; hold on; grid; xlim([1,N+1]);
plot(state_v(1,:),"r"); plot(state_v(2,:),"g"); plot(state_v(3,:),"b");
title("Velocity");legend("vx","vy","vz","Location","best");

figure; hold on; grid; xlim([1,N+1]);
plot(state_w(1,:),"r"); plot(state_w(2,:),"g"); plot(state_w(3,:),"b");
title("Angular Velocity");legend("wx","wy","wz","Location","best");

%% functions
function skew_mat = skew(vec)
    skew_mat = [  0     -vec(3)   vec(2);
                vec(3)     0     -vec(1);   
               -vec(2)   vec(1)     0  ];
end


function skew_3x12 = skew_4legs(leg_position)
     skew_3x12 = horzcat( skew(leg_position(1:3)), ...
                          skew(leg_position(4:6)), ...
                          skew(leg_position(7:9)), ...
                          skew(leg_position(10:12)) ...
                         );
end


function vec = vee(skew_mat)
    vec = [-skew_mat(2,3); skew_mat(1,3); -skew_mat(1,2)];
end


function Rot = exp_mat(skew_mat) 
    vec = vee(skew_mat);
    mag = sqrt( vec(1)^2 + vec(2)^2 + vec(3)^2 );
    idx = not(mag); % prevent division by zero
    Rot = eye(3) + sin(mag)/(mag+idx)*skew_mat + (1-cos(mag))/(mag^2+idx)*skew_mat^2;
end


function Rot = Exp_mat(vec)
    Rot = exp_mat(skew(vec));
end


function skew_mat = log_mat(Rot)
    if (trace(Rot) > 3.0-0.000001) % theta==0 (0.1deg resolution)
        skew_mat = 0.5*(Rot-Rot.'); %zeros(3,3);
    elseif (trace(Rot) < -1.0+0.000001) % theta==pi
        skew_vec = [pi*Rot(1,3)/sqrt(2*(1+Rot(3,3))), pi*Rot(2,3)/sqrt(2*(1+Rot(3,3))), pi*(1+Rot(3,3))/sqrt(2*(1+Rot(3,3)))];
        skew_mat = skew(skew_vec);
    else
        theta = acos( (trace(Rot)-1)/2 );
        skew_mat = theta/(2.0*sin(theta))*(Rot-Rot.');
    end
end


function vec = Log_mat(Rot)
    vec = vee( log_mat(Rot) );
end


function [rollout_p, rollout_v, rollout_R, rollout_w] = Euler_Integration(state_x,input_u)
    global N; global dt_mpc; global m_mpc; global IB IB_INV; global g; 
    global init_leg_position leg_position;
    global p0 v0 R0 w0;
    
    p = state_x(:,1);
    v = state_x(:,2);
    R = state_x(:,3:5);
    w = state_x(:,6);

    leg_position = init_leg_position - repmat(p-p0,4,1);
    F    = input_u(1:3) + input_u(4:6) + input_u(7:9) + input_u(10:12);
    tau  = cross(leg_position(1:3),input_u(1:3)) + cross(leg_position(4:6),input_u(4:6)) + cross(leg_position(7:9),input_u(7:9)) + cross(leg_position(10:12),input_u(10:12));
    v_dot = 1/m_mpc*F + g;
    w_dot = IB_INV * (R.'*tau - cross(w,IB*w));

    rollout_p = p + v*dt_mpc + 0.5*v_dot*dt_mpc^2;
    rollout_v = v + v_dot*dt_mpc;
    rollout_R = R*exp_mat(skew(w)*dt_mpc);
    rollout_w = w + w_dot*dt_mpc;
    
end


function objective = Objective_Function(state_p, state_v, state_R, state_w, input_U)
    global weight_p weight_v weight_psi weight_w weight_u weight_cone N;
    global p_ref v_ref R_ref w_ref;

    err_p = state_p(:,2:end) - p_ref;
    err_v = state_v(:,2:end) - v_ref;

    for k=1:N
        err_psi(:,k) = Log_mat( R_ref(:,3*k-2:3*k).' * state_R(:,3*(k+1)-2:3*(k+1)) );
        err_w(:,k)   = state_w(:,k+1) - Exp_mat(-err_psi(:,k))*w_ref(:,k);
    end

    obj_p   = sum(weight_p.*sum(err_p.*err_p,2)); 
    obj_v   = sum(weight_v.*sum(err_v.*err_v,2));
    obj_psi = sum(weight_psi.*sum(err_psi.*err_psi,2));
    obj_w   = sum(weight_w.*sum(err_w.*err_w,2));
    obj_U   = sum(weight_u.*sum(input_U.*input_U,2));

    objective = 0.5*(obj_p + obj_v + obj_psi + obj_w) + 0.5*obj_U;
end


function Jr = Right_Jacobian(err_psi)
    global N;
    
    I3 = eye(3);
    mag = norm(err_psi);
    idx = not(mag); % prevent division by zero
    X = skew( err_psi );
    Jr = I3 - (1-cos(mag))/(mag^2+idx)*X + (mag-sin(mag))/(mag^3+idx)*X^2;

end


function y = r_polar(apex, z)
    global mu_friction;
    y = 1/mu_friction * (apex - z);
end



function y = projection_cone_3x1(force_3x1, contact, LEG_NUM)
    global mu_friction; global fz_max; global fz_min; global isProjected;

    isProjected(LEG_NUM) = 1;
    y = zeros(3,1);

    if ( contact == 1 )
    
        radius = norm(force_3x1(1:2));
        idx = not(radius); % prevent division by zero (필요 없을듯? boolean으로 계산 시에는 필요.)
        r_polar_max = r_polar( (1+mu_friction^2)*fz_max , force_3x1(3) );
        r_polar_min = r_polar( (1+mu_friction^2)*fz_min , force_3x1(3) );
        
        section1 = (radius<mu_friction*fz_max) & (force_3x1(3)>fz_max);
        section2 = (radius>=r_polar_max) & (radius>=mu_friction*fz_max);
        section3 = (radius>mu_friction*force_3x1(3)) & (radius<r_polar_max) & (radius>r_polar_min);
        section4 = (radius>mu_friction*fz_min) & (radius<=r_polar_min);
        section5 = (radius<=mu_friction*fz_min) & (force_3x1(3)<fz_min);
      % section6 = (force_3x1(3)<=fz_max) & (force_3x1(3)>=fz_min) & (radius<=mu_friction*force_3x1(3));
        
        if section1
            y(3,1)   = fz_max;
            y(1:2,1) = force_3x1(1:2);
                
        elseif section2
            y(3,1)   = fz_max;
            y(1:2,1) = mu_friction*fz_max/radius * force_3x1(1:2);
            
        elseif section3
            y(3,1)   = (mu_friction*radius+force_3x1(3)) / (1+mu_friction^2);
            y(1:2,1) = mu_friction*y(3,1) * force_3x1(1:2)/radius;

        elseif section4
            y(3,1)   = fz_min;
            y(1:2,1) = mu_friction*fz_min/radius * force_3x1(1:2);
            
        elseif section5
            y(3,1)   = fz_min;
            y(1:2,1) = force_3x1(1:2);
            
        else
            isProjected(LEG_NUM) = 0;
            y = force_3x1;   
        end
    
    end
end



function y = projection_cone_12x1(force12x1, contact_4leg)
    y(1:3,1)   = projection_cone_3x1(force12x1(1:3,1), contact_4leg(1), 1);
    y(4:6,1)   = projection_cone_3x1(force12x1(4:6,1), contact_4leg(2), 2);
    y(7:9,1)   = projection_cone_3x1(force12x1(7:9,1), contact_4leg(3), 3);
    y(10:12,1) = projection_cone_3x1(force12x1(10:12,1), contact_4leg(4), 4);

end



function feedforward_constrained = Projected_GD(Quu, Quu_INV, Qu, input_u, contact_4leg)
    global isProjected; 
    
    % du의 projection operator는? Pc(u+du)-u
    del_u_init = projection_cone_12x1(input_u-Quu_INV*Qu, contact_4leg) - input_u;
    prev_isProjected = isProjected;

    du = del_u_init;
    du_prev = del_u_init;
    objective = 0.5*du.'*Quu*du + Qu.'*du;
    prev_objective = objective;
    initial_objective = objective;
    iter = 0;
    while 1
        iter = iter + 1;
        gradient = Quu*du + Qu;
        mag = norm(gradient);

        % backtracking
        Initial_Guess = 4000/mag;
        back_rate = 5;
        back_iter = 0;
        gamma = 0.001;
        while 1
            back_iter = back_iter + 1;
            L_smooth = back_rate^(back_iter-1) * 1/Initial_Guess;
            du_backtracking = projection_cone_12x1(input_u + (du-1/L_smooth*gradient), contact_4leg) - input_u;
            objective_backtracking = 0.5*du_backtracking.'*Quu*du_backtracking + Qu.'*du_backtracking;
            % back exit
            Gradient_Mapping = L_smooth * (du-du_backtracking);
            back_criteria = gamma/L_smooth * (Gradient_Mapping.'*Gradient_Mapping);
            back_condition1 = (objective - objective_backtracking) >= back_criteria;
            back_condition2 = back_iter >= 5;
            back_condition3 =  norm(Gradient_Mapping) < 1e-10;
            if (back_condition1 | back_condition2 | back_condition3) 
                break;
            end
        end
        
        % norm(Gradient_Mapping)
        if norm(Gradient_Mapping)<1e-10
            feedforward_constrained = du;
            isProjected = prev_isProjected;
            break
        end

        % next step update
        du = du_backtracking;
        objective = objective_backtracking;

        % main exit
        if iter == 1
            scale = abs(objective - prev_objective);
        end

        residual = ((objective - prev_objective)/scale)^2;
        exit_condition1 = residual <= 1e-3;
        exit_condition2 = iter>=10;
        if (exit_condition1 | exit_condition2) 
            break;
        end
        % update prev values 
        prev_objective = objective;
        du_prev = du;
        prev_isProjected = isProjected;
    end

    feedforward_constrained = du;
    
end



function free_space_Quu_INV = free_space_inverse(Quu,isProjected)
    
    % projection 된 곳을 지움!
    proj_expanded = ~ kron(isProjected,[1;1;1]);
    partitioned_Quu = proj_expanded*proj_expanded.'.*Quu;
    
    filled = sum(proj_expanded);
    non_zero_values = partitioned_Quu(partitioned_Quu ~= 0);
    compressed_Quu = reshape(non_zero_values,filled, filled);
    
    compressed_inv = inv(compressed_Quu);
    
    free_space_Quu_INV = zeros(12,12);
    
    num_i = 0; 
    for i = 1:4
        if isProjected(i)==0
            num_i = num_i+1;
        end
        num_j = 0;
        for j = 1:4
            if isProjected(j)==0
                num_j = num_j+1;
            end
    
            if (isProjected(i)==0) && (isProjected(j)==0)
                free_space_Quu_INV(3*i-2:3*i,3*j-2:3*j) = compressed_inv(3*num_i-2:3*num_i,3*num_j-2:3*num_j);
            end
    
        end
    end

end


function Feasibility_Check(input_U, contact_sequence)
    global isProjected cone_feasibility; global N;
    
    for k = 1:N
        projection_cone_12x1(input_U(:,k), contact_sequence(:,k));
        if sum( contact_sequence(:,k) .* isProjected ) > 0
            cone_feasibility = 0;
            break;         
        else
            cone_feasibility = 1;
        end
    end

end


function [del_u_feedforward, del_u_feedback_gain] = BackPass_DDP(input_U, contact_sequence) % for loop 안에서 k-1번째의state들을 인풋으로 받음.
    global weight_p weight_v weight_psi weight_w weight_u weight_cone;
    global terminal_weight_p terminal_weight_v terminal_weight_psi terminal_weight_w;
    global nom_error_p nom_error_v nom_error_psi nom_error_w nominal_p nominal_v nominal_R nominal_w;
    global p_ref v_ref R_ref w_ref;
    global p0 v0 R0 w0; 
    global state_p state_v state_R state_w;
    global N; global dt_mpc; global m_mpc; global IB IB_INV;
    global init_leg_position leg_position I3x12 I3; 
    global regularization_mu; global state_gap state_gap_init delta1 delta2;
    global Qu_data Quu_data Vx_data Vxx_data Qux_data;
    global isProjected cone_feasibility;

    % initialization
    fx = zeros(12,12); 
    fu = zeros(12,12);
    I12x12 = eye(12);

    terminal_weight = [terminal_weight_p; terminal_weight_v; terminal_weight_psi; terminal_weight_w];
    weight_x        = [weight_p; weight_v; weight_psi; weight_w];

    % state error calculations
    state_err_p = state_p(:,2:end) - p_ref;
    state_err_v = state_v(:,2:end) - v_ref;
    for k=1:N
        state_err_psi(:,k) = Log_mat( R_ref(:,3*k-2:3*k).' * state_R(:,3*(1+k)-2:3*(1+k)) );
        state_err_w(:,k)   = state_w(:,1+k) - Exp_mat(-state_err_psi(:,k))*w_ref(:,k);
    end
    init_err = [0;0;0];
    state_err_x = [init_err, state_err_p;
                   init_err, state_err_v;
                   init_err, state_err_psi;
                   init_err, state_err_w];

    Vxx = diag(terminal_weight) + regularization_mu*I12x12;
    Vx  = terminal_weight .* state_err_x(:,N+1) + Vxx * state_gap(:,N+1);

    for backward_step = 1 : N
        k = (N+1) - backward_step;
        leg_position = init_leg_position - repmat(nominal_p(:,k)-p0,4,1);

        fx(1:3,1:3) = I3; fx(1:3,4:6) = I3*dt_mpc; 
        fx(4:6,4:6) = I3;
        fx(7:9,7:9)   = Exp_mat(state_w(:,k)*dt_mpc).';
        fx(7:9,10:12) = Right_Jacobian(state_w(:,k)*dt_mpc)*dt_mpc;
        fx(10:12,7:9)   = IB_INV * skew(state_R(:,3*k-2:3*k).'*skew_4legs(leg_position)*input_U(:,k)) * dt_mpc;
        fx(10:12,10:12) = I3 - IB_INV * ( skew(state_w(:,k))*IB - skew(IB*state_w(:,k)) ) * dt_mpc;

        fu(1:3,:) = 0.5*dt_mpc^2/m_mpc * I3x12;
        fu(4:6,:) = dt_mpc/m_mpc * I3x12;
        fu(10:12,:) = IB_INV * state_R(:,3*k-2:3*k).' * skew_4legs(leg_position) * dt_mpc;

        fu = fu.*kron(contact_sequence(:,k).',[1 1 1]);
        
        Qx = weight_x.*state_err_x(:,k) + fx.'*Vx; % weight_x.'*로 해버려서 스칼라로 더해져버려서 이상하게 나왔었음.. C++로 할땐 에러가 뜨겠지만 쨌든 diag(weight_x)의 의미가 맞다.
        Qu = weight_u.*input_U(:,k) + fu.'*Vx;
        Qxx = diag(weight_x) + fx.'*Vxx*fx;
        Quu = diag(weight_u) + fu.'*Vxx*fu + regularization_mu*I12x12; Quu_INV = inv(Quu);
        Qux = fu.'*Vxx*fx;
        
        if cone_feasibility == 0
            del_u_feedforward(:,k) = -Quu_INV*Qu; 
            del_u_feedback_gain(:,12*k-11:12*k) = -Quu_INV*Qux; 
        else
            del_u_feedforward(:,k) = Projected_GD(Quu, Quu_INV, Qu, input_U(:,k), contact_sequence(:,k));
            del_u_feedback_gain(:,12*k-11:12*k) = -free_space_inverse(Quu,isProjected)*Qux;
        end

        Vxx = Qxx + Qux.' * del_u_feedback_gain(:,12*k-11:12*k) + regularization_mu*I12x12; 
        Vx  = Qx  + Qux.' * del_u_feedforward(:,k) + Vxx*state_gap(:,k);  
        
        Qu_data(:,k) = Qu; Vx_data(:,k) = Vx; Qux_data(:,12*k-11:12*k) = Qux;
        Quu_data(:,12*k-11:12*k) = Quu; Vxx_data(:,12*k-11:12*k) = Vxx;
    end
    
end


function [optimal_U, optimal_p, optimal_v, optimal_R, optimal_w] = ForwardPass_DDP(del_u_feedforward, del_u_feedback_gain, input_U, contact_sequence)
    global N; global dt_mpc; global m_mpc; global IB IB_INV; global g; 
    global init_leg_position leg_position;
    global p_ref v_ref R_ref w_ref;
    global p0 v0 R0 w0;
    global state_p state_v state_R state_w;
    global backtracking_alpha step_size_alpha;
    global nom_error_p nom_error_v nom_error_psi nom_error_w nominal_p nominal_v nominal_R nominal_w;
    global state_gap state_gap_init delta1 delta2;
    global Qu_data Quu_data Vx_data Vxx_data Qux_data;
    global isProjected;
    
    % alpha step to the input direction
    del_x = zeros(12,1);

    optimal_u = input_U(:,1) ...
              + step_size_alpha * del_u_feedforward(:,1) ...
              + del_u_feedback_gain(:,1:12)*del_x;
    
    optimal_U(:,1) = projection_cone_12x1(optimal_u,contact_sequence(:,1));

    % alpha step to the state direction
    [rollout_p, rollout_v, rollout_R, rollout_w] = Euler_Integration([state_p(:,1),state_v(:,1),state_R(:,3*1-2:3*1),state_w(:,1)], optimal_U(:,1));
    
    optimal_p(:,1)   = rollout_p + (step_size_alpha-1)*state_gap(1:3,1+1);
    optimal_v(:,1)   = rollout_v + (step_size_alpha-1)*state_gap(4:6,1+1);
    optimal_R(:,1:3) = rollout_R * Exp_mat((step_size_alpha-1)*state_gap(7:9,1+1));
    optimal_w(:,1)   = rollout_w + (step_size_alpha-1) * state_gap(10:12,1+1);

    % expected improvement
    delta1 = 0; delta2 = 0;
    delta1 = delta1 + Qu_data(:,1).' * del_u_feedforward(:,1) ...
                    + state_gap(:,1).'*(Vx_data(:,1) - Vxx_data(:,1:12)*del_x);
    delta2 = delta2 + del_u_feedforward(:,1).'*Quu_data(:,1:12)*del_u_feedforward(:,1) ...
                    + state_gap(:,1).'*(2*Vxx_data(:,1:12)*del_x - Vxx_data(:,1:12)*state_gap(:,1));

    for k = 2:N
        % alpha step to the input direction
        del_x = [optimal_p(:,k-1) - state_p(:,k);
                 optimal_v(:,k-1) - state_v(:,k);
                 Log_mat(state_R(:,3*(k)-2:3*(k)).'*optimal_R(:,3*(k-1)-2:3*(k-1))); 
                 optimal_w(:,k-1) - (state_R(:,3*(k)-2:3*(k)).'*optimal_R(:,3*(k-1)-2:3*(k-1))).'*state_w(:,k)];


        optimal_u = input_U(:,k) ...
                  + step_size_alpha * del_u_feedforward(:,k) ...
                  + del_u_feedback_gain(:,12*k-11:12*k)*del_x;

        optimal_U(:,k) = projection_cone_12x1(optimal_u,contact_sequence(:,k));
        
        % alpha step to the state direction
        [rollout_p, rollout_v, rollout_R, rollout_w] = Euler_Integration([optimal_p(:,k-1),optimal_v(:,k-1),optimal_R(:,3*(k-1)-2:3*(k-1)),optimal_w(:,k-1)], optimal_U(:,k));
        
        optimal_p(:,k)         = rollout_p + (step_size_alpha-1)*state_gap(1:3,1+k);
        optimal_v(:,k)         = rollout_v + (step_size_alpha-1)*state_gap(4:6,1+k);
        optimal_R(:,3*k-2:3*k) = rollout_R * Exp_mat((step_size_alpha-1)*state_gap(7:9,1+k));
        optimal_w(:,k)         = rollout_w + (step_size_alpha-1)*state_gap(10:12,1+k);

        % expected improvement 여기의 state_gap은 index +1 안하는게 맞다. 변수들 다 0~N에 맞는 인덱스로 되어 있음.
        delta1 = delta1 + Qu_data(:,k).' * del_u_feedforward(:,k) ...
                        + state_gap(:,k).'*(Vx_data(:,k) - Vxx_data(:,12*k-11:12*k)*del_x);
        delta2 = delta2 + del_u_feedforward(:,k).'*Quu_data(:,12*k-11:12*k)*del_u_feedforward(:,k) ...
                        + state_gap(:,k).'*(2*Vxx_data(:,12*k-11:12*k)*del_x - Vxx_data(:,12*k-11:12*k)*state_gap(:,k));

    end

end


%%
    % % Forward Pass
    % sub_iter = 0;
    % step_size_alpha = 0.8;  
    % while 1
    %     sub_iter = sub_iter + 1;
    %     state_gap = state_gap_init;
    %     % step_size_alpha = 0.5*step_size_alpha;
    %     % step_size_alpha = 1.2^(-(sub_iter-1));
    % 
    %     [optimal_U, optimal_p, optimal_v, optimal_R, optimal_w] = ForwardPass_DDP(del_u_feedforward, del_u_feedback_gain, input_U);
    %     objective_forwardpass = Objective_Function([p0,optimal_p], [v0,optimal_v], [R0, optimal_R], [w0,optimal_w], optimal_U);
    % 
    %     expected_improvement = delta1*step_size_alpha + 0.5*delta2*step_size_alpha^2;
    %     delta1
    %     delta2
    %     % step_size_alpha = 0.5*(step_size_alpha + max(min(-delta1/delta2,1.0),0) )
    % 
    %     if expected_improvement <= 0
    %         b = 0;
    %     else
    %         b = 2e-2;
    %     end
    % 
    %     % step_size_alpha = exp(-(1/2e5* 0.65 *expected_improvement )^2) %1/(1+exp(1/1e5*expected_improvement));
    %     back_condition1 = objective - prev_objective < b*expected_improvement; 
    %     back_condition2 = sub_iter >= MAX_SUB_ITER;
    % 
    %     if (back_condition1 | back_condition2) 
    %         fprintf("sub_iter: "+string(sub_iter)+"\n");
    %         fprintf("expected_improvement: "+string(expected_improvement)+"\n");
    %         break;
    %     end
    % end

    %%




