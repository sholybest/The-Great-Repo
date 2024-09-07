%%
clear; clc;

%%
global N; global dt_mpc; global m_mpc; global IB IB_INV; 
global init_leg_position leg_position I3x12 I3; global g; 
global mu_friction; global fz_max; global fz_min;
global weight_p weight_v weight_psi weight_w weight_u weight_cone;
global terminal_weight_p terminal_weight_v terminal_weight_psi terminal_weight_w;
global backtracking_alpha;
global p_ref v_ref R_ref w_ref;
global nom_error_p nom_error_v nom_error_psi nom_error_w nominal_p nominal_v nominal_R nominal_w;
global p0 v0 R0 w0;
global cone_mat dd_cone weight_linear_cone;
global regularization_mu;
global Qu_data Quu_data Vx_data Vxx_data Qux_data;
global isProjected; 

%% Parameters
N = 25;
dt_mpc = 0.02;
m_mpc = 40;
IB = diag([0.4 2.1 2.1]);
IB_INV = inv(IB);
mu_friction = 0.6; fz_max = 666; fz_min = 10;
I3 = eye(3);

isProjected = [0;0;0;0];
regularization_mu = 0;
backtracking_alpha = 1.0;
back_ratio = 2.5; %2.5
MAX_ITER = 50;
MAX_BACK_ITER = 10; % 6,10
epsilon_main = 1e-6;

% weight_p   = [10;10;5000]; % 1;1;50?
% weight_v   = [10;10;10];  
% weight_psi = [200;200;200]; % 5;5;5 or 60;60;60
% weight_w   = [1e1;1e1;1e1];
% weight_u   = [1;1;1; 1;1;1; 1;1;1; 1;1;1]*1e-3;

weight_p   = [1;1;50]; % 1;1;50?
weight_v   = [1;1;0.1];  
weight_psi = [1;1;1]; % 5;5;5 or 60;60;60
weight_w   = [1e-2;1e-2;1e-2];
weight_u   = [1;1;1; 1;1;1; 1;1;1; 1;1;1]*1e-6;

terminal_weight_p   = [1;1;500]; %10.0*weight_p;
terminal_weight_v   = 1.0*weight_v;
terminal_weight_psi = [1;1;1];
terminal_weight_w   = 1.0*weight_w;

%% Initialization
% r1 = [ 0.317;  0.12; -0.48];
% r2 = [ 0.317; -0.12; -0.48];
% r3 = [-0.317;  0.12; -0.48];
% r4 = [-0.317; -0.12; -0.48];

r1 = [ 0.;  0.12; -0.48];
r2 = [ 0.; -0.12; -0.48];
r3 = [-0.;  0.12; -0.48];
r4 = [-0.; -0.12; -0.48];

init_leg_position = [r1; r2; r3; r4]; % 발 위치가 변할 경우 옆으로 이어붙이기.

I3x12 = repmat(eye(3),1,4);
g = [0;0;-9.81];

p_ref = repmat([0;0;0.5],1,N);
v_ref = repmat([0;0;0],1,N);
R_ref = repmat(eye(3),1,N); 
w_ref = repmat([0;0;0],1,N);

p0 = [0;0;0.48];
v0 = [0.01;-0.01;-0.01];
R0 = rotz(.1)*roty(-.8)*rotx(.1);
w0 = [0.01;-0.01;0.01];

% p0 = [0;0.1;0.27];
% v0 = [0.05;-0.05;-0.04];
% R0 = rotz(12)*roty(15)*rotx(14);
% w0 = [0.1;0.1;0.1];

state_p = zeros(3,N);
state_v = zeros(3,N);
state_R = zeros(3,3*N); % 3x3N
state_w = zeros(3,N);

error_p   = zeros(3,N);
error_v   = zeros(3,N);
error_psi = zeros(3,N);
error_w   = zeros(3,N);

%% DDP iterations
% initial nominal trajectory
% initial_input_U = repmat([0;0;100; 0;0;100; 0;0;100; 0;0;100],1,N); 
initial_input_U = repmat([0;0;100; 0;0;0; 0;0;0; 0;0;100],1,N);
% contact_sequence = repmat([1;1;1;1],1,N);
contact_sequence = repmat([1;0;0;1],1,N);

% load input_U.mat;
% initial_input_U = input_U;

[nom_error_p, nom_error_v, nom_error_psi, nom_error_w, nominal_p, nominal_v, nominal_R, nominal_w] = Model_Update(initial_input_U);
objective = Objective_Function(nom_error_p, nom_error_v, nom_error_psi, nom_error_w, initial_input_U); obj_data(1) = objective;

tic;

prev_objective = 0;
initial_objective = objective;
minimum_objective = objective;
input_U = initial_input_U;
iter = 0;
while 1
    iter = iter + 1; fprintf("\n\niteration "+string(iter)+":\n");

    % Backward pass
    [del_u_feedforward, del_u_feedback_gain] = BackPass_DDP(input_U, contact_sequence);

    % Forward Pass
    back_iter = 0;
    while 1
        back_iter = back_iter + 1;
        backtracking_alpha = back_ratio^(-(back_iter-1));
        [optimal_U, opt_error_p, opt_error_v, opt_error_psi, opt_error_w, optimal_p, optimal_v, optimal_R, optimal_w] = ForwardPass_DDP(del_u_feedforward, del_u_feedback_gain, input_U, contact_sequence);
        objective_forwardpass = Objective_Function(opt_error_p, opt_error_v, opt_error_psi, opt_error_w, optimal_U);
        
        % forwardpass 함수 내의 마지막 스텝 del_x도 저장하는거 잊지 말기! 
        Qux = Qux_data(:,12*N-11:12*N);
        Quu = Quu_data(:,12*N-11:12*N);
        Qu = Qu_data(:,N);
        save Qux.mat Qux;
        save Quu.mat Quu;
        save Qu.mat Qu;

        back_condition1 = objective > objective_forwardpass;
        back_condition2 = back_iter >= MAX_BACK_ITER;
        if (back_condition1 | back_condition2) 
            fprintf("back_iter: "+string(back_iter)+"\n");
            break;
        end
    end

    % Next step update
    input_U = optimal_U;
    Update_Nominal_States(optimal_p,optimal_v,optimal_R,optimal_w, opt_error_p,opt_error_v,opt_error_psi,opt_error_w);
    prev_objective = objective;
    objective = objective_forwardpass; obj_data(iter+1) = objective;

    % save minimal value
    if (objective <= minimum_objective)
        minimum_objective = objective;
        minimum_U = input_U;
    end

    % exit conditions
    residual = ((prev_objective-objective)/initial_objective)^2;
    fprintf("residual: "+string(residual)+"\n");
    exit_condition1 = residual <= epsilon_main;
    exit_condition2 = 0;%~ back_condition1; % backtracking failed
    exit_condition3 = iter >= MAX_ITER;
    
    if (exit_condition1 | exit_condition2 | exit_condition3)
        fprintf("exit_condition1: "+string(double(exit_condition1))+"\n");
        fprintf("exit_condition2: "+string(double(exit_condition2))+"\n");
        fprintf("exit_condition3: "+string(double(exit_condition3))+"\n");
        break;
    end
   
end

toc;

%% plots
figure; plot(0:iter,obj_data); grid;
title("Objective Function"); xlabel("iterations"); ylabel("Objective");

figure; plot(vecnorm(opt_error_psi*180/pi)); grid; 
xlim([1 N]); %ylim([0 25]);
title("Rotational error"); xlabel("steps"); ylabel("rotational error (deg)");

figure; hold on; grid; xlim([1,N+1]);
plot(optimal_p(1,:),"r"); plot(optimal_p(2,:),"g"); plot(optimal_p(3,:),"b");
title("Position"); legend("px","py","pz","Location","northwest");

figure; hold on; grid; xlim([1,N+1]);
plot(optimal_v(1,:),"r"); plot(optimal_v(2,:),"g"); plot(optimal_v(3,:),"b");
title("Velocity");legend("vx","vy","vz","Location","best");

figure; hold on; grid; xlim([1,N+1]);
plot(optimal_w(1,:),"r"); plot(optimal_w(2,:),"g"); plot(optimal_w(3,:),"b");
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


function [err_p, err_v, err_psi, err_w, update_p, update_v, update_R, update_w] = Model_Update(input_U)
    global N; global dt_mpc; global m_mpc; global IB IB_INV; global g; 
    global init_leg_position leg_position;
    global p_ref v_ref R_ref w_ref;
    global p0 v0 R0 w0;

    % First update (obtain x1's from x0's)
    leg_position = init_leg_position;
    F0    = input_U(1:3,1) + input_U(4:6,1) + input_U(7:9,1) + input_U(10:12,1);
    tau0  = cross(leg_position(1:3),input_U(1:3,1)) + cross(leg_position(4:6),input_U(4:6,1)) + cross(leg_position(7:9),input_U(7:9,1)) + cross(leg_position(10:12),input_U(10:12,1));
    v0_dot = 1/m_mpc*F0 + g;
    w0_dot = IB_INV * (R0.'*tau0 - cross(w0,IB*w0));

    update_p = zeros(3,N); % setzero?
    update_v = zeros(3,N);
    update_R = zeros(3,3*N);
    update_w = zeros(3,N);

    update_p(:,1) = p0 + v0*dt_mpc + 0.5*v0_dot*dt_mpc^2;
    update_v(:,1) = v0 + v0_dot*dt_mpc;
    update_R(:,1:3) = R0*exp_mat(skew(w0)*dt_mpc);
    update_w(:,1) = w0 + w0_dot*dt_mpc;

    err_psi(:,1) = Log_mat(R_ref(:,1:1+2).'*update_R(:,1:3));
    err_w(:,1)   = update_w(:,1) - Exp_mat(-err_psi(:,1))*w_ref(:,1);

    % Model updates (obtain x2 ~ xN)
    for k = 2:N
        leg_position = init_leg_position - repmat(update_p(:,k-1)-p0,4,1);
        F_sum  = input_U(1:3,k) + input_U(4:6,k) + input_U(7:9,k) + input_U(10:12,k);
        tau    = cross(leg_position(1:3),input_U(1:3,k)) + cross(leg_position(4:6),input_U(4:6,k)) + cross(leg_position(7:9),input_U(7:9,k)) + cross(leg_position(10:12),input_U(10:12,k));
        v_dot = 1/m_mpc*F_sum + g;
        w_dot = IB_INV * (update_R(:,3*(k-1)-2:3*(k-1)).'*tau - cross(update_w(:,k-1),IB*update_w(:,k-1)));
        
        update_p(:,k)   = update_p(:,k-1) + update_v(:,k-1)*dt_mpc + 0.5*v_dot*dt_mpc^2;
        update_v(:,k)   = update_v(:,k-1) + v_dot*dt_mpc;
        update_R(:,3*k-2:3*k) = update_R(:,3*(k-1)-2:3*(k-1))*exp_mat(skew(update_w(:,k-1))*dt_mpc); 
        update_w(:,k) = update_w(:,k-1) + w_dot*dt_mpc;
        
        err_psi(:,k) = Log_mat( R_ref(:,3*k-2:3*k).' * update_R(:,3*k-2:3*k) );
        err_w(:,k)   = update_w(:,k) - Exp_mat(-err_psi(:,k))*w_ref(:,k);
    end
 
    err_p = update_p - p_ref;
    err_v = update_v - v_ref;
end

function first_column = COL1(mat)
    first_column = mat(:,1);
end

function second_column = COL2(mat)
    second_column = mat(:,2);
end

function third_column = COL3(mat)
    third_column = mat(:,3);
end


function objective = Objective_Function(err_p, err_v, err_psi, err_w, input_U)
    global weight_p weight_v weight_psi weight_w weight_u weight_cone N;
    global cone_mat dd_cone weight_linear_cone;

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


function Jr_inv = Right_Jacobian_INV(err_psi)
    global N;
    
    I3 = eye(3);
    mag = norm(err_psi);
    idx = not(mag); % prevent division by zero
    X = skew( err_psi );
    Jr_inv = I3 + 0.5*X + (1/(mag^2+idx)-(1+cos(mag))/(2*mag*sin(mag)+idx))*X^2;

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
        Initial_Guess = 2000/mag;
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
        exit_condition2 = iter>=5;
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



function free_space_Quu_INV = free_space_inverse(Quu)
    global isProjected;
    
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



function [del_u_feedforward, del_u_feedback_gain] = BackPass_DDP(input_U, contact_sequence) % for loop 안에서 k-1번째의state들을 인풋으로 받음.
    global weight_p weight_v weight_psi weight_w weight_u weight_cone;
    global terminal_weight_p terminal_weight_v terminal_weight_psi terminal_weight_w;
    global nom_error_p nom_error_v nom_error_psi nom_error_w nominal_p nominal_v nominal_R nominal_w;
    global p0 v0 R0 w0; 
    global N; global dt_mpc; global m_mpc; global IB IB_INV;
    global init_leg_position leg_position I3x12 I3; 
    global cone_mat dd_cone weight_linear_cone;
    global regularization_mu;
    global Qu_data Quu_data Vx_data Vxx_data Qux_data;
    global isProjected; 

    % nominal trajectories
    init_err = [0;0;0];
    nominal_err_x = [init_err, nom_error_p;
                     init_err, nom_error_v;
                     init_err, nom_error_psi;
                     init_err, nom_error_w];
    nominal_state_R = [R0, nominal_R];
    nominal_state_w = [w0, nominal_w];
    
    fx = zeros(12,12); 
    fu = zeros(12,12);
    I12x12 = eye(12);
    diff_cone_error = zeros(12,12);

    terminal_weight = [terminal_weight_p; terminal_weight_v; terminal_weight_psi; terminal_weight_w];
    weight_x        = [weight_p; weight_v; weight_psi; weight_w];

    % for terminal cost (xN)
    Vx  = terminal_weight .* nominal_err_x(:,N+1);
    Vxx = diag(terminal_weight);

    for backward_step = 1 : N
        k = (N+1) - backward_step;
        leg_position = init_leg_position - repmat(nominal_p(:,k)-p0,4,1);

        fx(1:3,1:3) = I3; fx(1:3,4:6) = I3*dt_mpc; 
        fx(4:6,4:6) = I3;
        fx(7:9,7:9)   = Exp_mat(nominal_state_w(:,k)*dt_mpc).';
        fx(7:9,10:12) = Right_Jacobian(nominal_state_w(:,k)*dt_mpc)*dt_mpc;
        fx(10:12,7:9)   = IB_INV * skew(nominal_state_R(:,3*k-2:3*k).'*skew_4legs(leg_position)*input_U(:,k)) * dt_mpc;
        fx(10:12,10:12) = I3 - IB_INV * ( skew(nominal_state_w(:,k))*IB - skew(IB*nominal_state_w(:,k)) ) * dt_mpc;

        fu(1:3,:) = 0.5*dt_mpc^2/m_mpc * I3x12;
        fu(4:6,:) = dt_mpc/m_mpc * I3x12;
        fu(10:12,:) = IB_INV * nominal_state_R(:,3*k-2:3*k).' * skew_4legs(leg_position) * dt_mpc;

        fu = fu.*kron(contact_sequence(:,k).',[1 1 1]);
        
        Qx = weight_x.*nominal_err_x(:,k) + fx.'*Vx; 
        Qu = weight_u.*input_U(:,k) + fu.'*Vx;
        Qxx = diag(weight_x) + fx.'*Vxx*fx;
        Quu = diag(weight_u) + fu.'*Vxx*fu + regularization_mu*I12x12; Quu_INV = inv(Quu);
        Qux = fu.'*Vxx*fx;        

        % del_u_feedforward(:,k) = -Quu_INV*Qu; 
        % del_u_feedback_gain(:,12*k-11:12*k) = -Quu_INV*Qux;

        del_u_feedforward(:,k) = Projected_GD(Quu, Quu_INV, Qu, input_U(:,k), contact_sequence(:,k));
        del_u_feedback_gain(:,12*k-11:12*k) = -free_space_inverse(Quu)*Qux;

        Vx  = Qx  + Qux.' * del_u_feedforward(:,k);
        Vxx = Qxx + Qux.' * del_u_feedback_gain(:,12*k-11:12*k);
        
        Qu_data(:,k) = Qu; Vx_data(:,k) = Vx; Qux_data(:,12*k-11:12*k) = Qux;
        Quu_data(:,12*k-11:12*k) = Quu; Vxx_data(:,12*k-11:12*k) = Vxx;
    end
    
end



function [optimal_U, opt_error_p, opt_error_v, opt_error_psi, opt_error_w, optimal_p, optimal_v, optimal_R, optimal_w] ...
                    = ForwardPass_DDP(del_u_feedforward, del_u_feedback_gain, input_U, contact_sequence)
    global N; global dt_mpc; global m_mpc; global IB IB_INV; global g; 
    global init_leg_position leg_position;
    global p_ref v_ref R_ref w_ref;
    global p0 v0 R0 w0;
    global backtracking_alpha;
    global nom_error_p nom_error_v nom_error_psi nom_error_w nominal_p nominal_v nominal_R nominal_w;
    global Qu_data Quu_data Vx_data Vxx_data Qux_data;
    global isProjected; 

    del_x = zeros(12,1);
    optimal_U(:,1) = input_U(:,1) ...
                    + backtracking_alpha * del_u_feedforward(:,1) ...
                    + del_u_feedback_gain(:,1:12)*del_x;

    % First update (obtain x1's from x0's)
    leg_position = init_leg_position;
    F0    = optimal_U(1:3,1) + optimal_U(4:6,1) + optimal_U(7:9,1) + optimal_U(10:12,1);
    tau0  = cross(leg_position(1:3),optimal_U(1:3,1)) + cross(leg_position(4:6),optimal_U(4:6,1)) + cross(leg_position(7:9),optimal_U(7:9,1)) + cross(leg_position(10:12),optimal_U(10:12,1));
    v0_dot = 1/m_mpc*F0 + g;
    w0_dot = IB_INV * (R0.'*tau0 - cross(w0,IB*w0));

    optimal_p(:,1) = p0 + v0*dt_mpc + 0.5*v0_dot*dt_mpc^2;
    optimal_v(:,1) = v0 + v0_dot*dt_mpc;
    optimal_R(:,1:3) = R0*exp_mat(skew(w0)*dt_mpc);
    optimal_w(:,1) = w0 + w0_dot*dt_mpc;

    opt_error_psi(:,1) = Log_mat(R_ref(:,1:1+2).'*optimal_R(:,1:3));
    opt_error_w(:,1)   = optimal_w(:,1) - Exp_mat(-opt_error_psi(:,1))*w_ref(:,1);

    % Model updates (obtain x2 ~ xN)
    for k = 2:N
        del_x = [optimal_p(:,k-1) - nominal_p(:,k-1);
                 optimal_v(:,k-1) - nominal_v(:,k-1);
                 Log_mat(nominal_R(:,3*(k-1)-2:3*(k-1)).'*optimal_R(:,3*(k-1)-2:3*(k-1))); %opt_error_psi(:,k-1) - nom_error_psi(:,k-1);
                 optimal_w(:,k-1) - (optimal_R(:,3*(k-1)-2:3*(k-1)).'*nominal_R(:,3*(k-1)-2:3*(k-1))).'*nominal_w(:,k-1)];
                

        optimal_u = input_U(:,k) ...
                + backtracking_alpha * del_u_feedforward(:,k) ...
                + del_u_feedback_gain(:,12*k-11:12*k)*del_x;

        optimal_U(:,k) = projection_cone_12x1(optimal_u,contact_sequence(:,k));

        leg_position = init_leg_position - repmat(optimal_p(:,k-1)-p0,4,1);
        F_sum  = optimal_U(1:3,k) + optimal_U(4:6,k) + optimal_U(7:9,k) + optimal_U(10:12,k);
        tau    = cross(leg_position(1:3),optimal_U(1:3,k)) + cross(leg_position(4:6),optimal_U(4:6,k)) + cross(leg_position(7:9),optimal_U(7:9,k)) + cross(leg_position(10:12),optimal_U(10:12,k));
        v_dot = 1/m_mpc*F_sum + g;
        w_dot = IB_INV * (optimal_R(:,3*(k-1)-2:3*(k-1)).'*tau - cross(optimal_w(:,k-1),IB*optimal_w(:,k-1)));
        
        optimal_p(:,k)   = optimal_p(:,k-1) + optimal_v(:,k-1)*dt_mpc + 0.5*v_dot*dt_mpc^2;
        optimal_v(:,k)   = optimal_v(:,k-1) + v_dot*dt_mpc;
        optimal_R(:,3*k-2:3*k) = optimal_R(:,3*(k-1)-2:3*(k-1))*exp_mat(skew(optimal_w(:,k-1))*dt_mpc); 
        optimal_w(:,k) = optimal_w(:,k-1) + w_dot*dt_mpc;
        
        opt_error_psi(:,k) = Log_mat( R_ref(:,3*k-2:3*k).' * optimal_R(:,3*k-2:3*k) );
        opt_error_w(:,k)   = optimal_w(:,k) - Exp_mat(-opt_error_psi(:,k))*w_ref(:,k);
    end
 
    opt_error_p = optimal_p - p_ref;
    opt_error_v = optimal_v - v_ref;
    save del_x.mat del_x;
end

function Update_Nominal_States(optimal_p, optimal_v, optimal_R, optimal_w, opt_error_p, opt_error_v, opt_error_psi, opt_error_w)
global nom_error_p nom_error_v nom_error_psi nom_error_w nominal_p nominal_v nominal_R nominal_w;

    nominal_p = optimal_p;
    nominal_v = optimal_v;
    nominal_R = optimal_R;
    nominal_w = optimal_w;
    nom_error_p = opt_error_p;
    nom_error_v = opt_error_v;
    nom_error_psi = opt_error_psi;
    nom_error_w = opt_error_w;

end



