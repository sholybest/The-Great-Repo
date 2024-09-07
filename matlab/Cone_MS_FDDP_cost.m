%%
clear; clc;

%%
global N; global dt_mpc; global m_mpc; global IB IB_INV; 
global init_leg_position I3x12 I3 I12; global g; 
global mu_friction; global fz_max; global fz_min;
global weight_p weight_v weight_psi weight_w weight_u weight_cone;
global terminal_weight_p terminal_weight_v terminal_weight_psi terminal_weight_w;
global step_size_alpha;
global p_ref v_ref R_ref w_ref;
global nominal_p nominal_v nominal_R nominal_w;
global p0 v0 R0 w0;
global state_p state_v state_R state_w;
global regularization_mu;
global state_gap state_gap_init delta1 delta2;
global Vx_data Vxx_data;
global weight_state_bound x_lower x_upper;

%% Parameters
N = 25;
dt_mpc = 0.02;
m_mpc = 40;
IB = diag([0.4 2.1 2.1]);
IB_INV = inv(IB);
mu_friction = 0.6; fz_max = 666; fz_min = 10;
I3 = eye(3); I12 = eye(12);

Qu_data = zeros(12,N);
Quu_data = zeros(12,12*N);
Vx_data = zeros(12,N);
Vxx_data = zeros(12,12*N);
state_gap = zeros(12,N+1);
delta1 = 0; delta2 = 0;
step_size_alpha = 0.5;
regularization_mu = 1e-6;

MAX_ITER = 5000;
epsilon_gap = 1e-6;

% weight_p   = [50;50;50];
% weight_v   = [1;1;1];  
% weight_psi = [1;1;1]; 
% weight_w   = [1e-2;1e-2;1e-2];
% weight_u   = [1;1;1; 1;1;1; 1;1;1; 1;1;1]*1e-4;
% 
% terminal_weight_p   = 1*weight_p;
% terminal_weight_v   = 1*weight_v;
% terminal_weight_psi = 1*weight_psi;
% terminal_weight_w   = 1*weight_w;
% 
% weight_cone = ones(12,1)*10;
% 
% weight_state_bound = [0;0;0;
%                       0;0;1e3;
%                       0;0;0;
%                       10;10;10];

weight_p   = [50;50;50];
weight_v   = [1;1;1];  
weight_psi = [20;20;20]; 
weight_w   = [1e-2;1e-2;1e-2];
weight_u   = [1;1;1; 1;1;1; 1;1;1; 1;1;1]*0e-6;

terminal_weight_p   = 1*weight_p;
terminal_weight_v   = 1*weight_v;
terminal_weight_psi = 1*weight_psi;
terminal_weight_w   = 1*weight_w;

weight_cone = ones(12,1)*1e6;

% weight_state_bound = [0;0;0;
%                       0;0;1e3;
%                       0;0;0;
%                       10;10;10];

weight_state_bound = [0;0;0;
                      0;0;0;
                      0;0;0;
                      0;0;0];

x_lower = [0; 0; 0;   % position
           0; 0; -0.1;   % velocity
           0; 0; 0;  % rotational error
           -0.5; -0.5; -0.5]; % angular veolcity 

x_upper = [0; 0; 0;   % position
           0; 0; 0.1;   % velocity
           0; 0; 0;  % rotational error
           0.5; 0.5; 0.5]; % angular veolcity 

%% Initialization
H0 = 0.44;
r1 = [ 0.317;  0.12; -H0];
r2 = [ 0.317; -0.12; -H0];
r3 = [-0.317;  0.12; -H0];
r4 = [-0.317; -0.12; -H0];

% r1 = [ 0.;  0.12; -H0];
% r2 = [ 0.; -0.12; -H0];
% r3 = [-0.;  0.12; -H0];
% r4 = [-0.; -0.12; -H0];

init_leg_position = [r1; r2; r3; r4]; % maxCoeff발 위치가 변할 경우 옆으로 이어붙이기.

I3x12 = repmat(eye(3),1,4);
g = [0;0;-9.81];

p_ref = repmat([0.05;0;0.5],1,N);
v_ref = repmat([0;0;0],1,N);
R_ref = repmat(eye(3,3),1,N); 
w_ref = repmat([0;0;0],1,N);

% p0 = [0;0;H0];
% v0 = [0.01;-0.1;-0.01];
% R0 = rotz(0.1)*roty(-2)*rotx(1);
% w0 = [0.01;-0.01;0.01];

% p0 = [0.0;0.0;H0];
% v0 = [-0.189;0.0709;-0.121];
% R0 = rotz(0.09*180/pi)*roty(0.03*180/pi)*rotx(0.0006*180/pi);
% w0 = [5.9;0.46;0.033];

p0 = [0;-0.01;H0];
v0 = [0.05;-0.01;-0.04];
R0 = rotz(15)*roty(12)*rotx(14);
w0 = [0.1;0.1;0.1];

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

% initial nominal trajectory
initial_input_U = repmat([0;0;100; 0;0;100; 0;0;100; 0;0;100],1,N); 
% initial_input_U = repmat([0;0;200; 0;0;0; 0;0;0; 0;0;200],1,N);
contact_sequence = repmat([1;1;1;1],1,N);
% contact_sequence = repmat([1;1;0;0],1,N);

%% quasi static input
% w1 = 1;
% w2 = 1;
% 
% A1 = [I3 I3 I3 I3];
% A2 = IB_INV*R0.'*[skew(r1) skew(r2) skew(r3) skew(r4)];
% 
% b1 = -m_mpc*(v0/dt_mpc+g);
% b2 = IB_INV*skew(w0)*IB*w0 - w0/dt_mpc;
% % b1 = -m_mpc*g;
% % b2 = IB_INV*skew(w0)*IB*w0;
% 
% contact0 = contact_sequence(:,1);
% C_ = [1  0 -mu_friction;
%     -1  0 -mu_friction;
%      0  1 -mu_friction;
%      0 -1 -mu_friction;
%      0  0 -1;
%      0  0  1];
% d_ = [0;0;0;0;-fz_min; fz_max];
% 
% A = blkdiag(C_,C_,C_,C_);
% b = [contact0(1)*d_;contact0(2)*d_;contact0(3)*d_;contact0(4)*d_];
% 
% H = w1*A1.'*A1 + w2*A2.'*A2;
% f = -(w1*A1.'*b1 + w2*A2.'*b2);
% 
% u = quadprog(H,f,A,b);

%% DDP iterations
nominal_p = repmat(p0,1,N);
nominal_v = repmat(v0,1,N);
nominal_R = repmat(R0,1,N);
nominal_w = repmat(w0,1,N);

state_p = [p0, nominal_p];
state_v = [v0, nominal_v];
state_R = [R0, nominal_R];
state_w = [w0, nominal_w];

% state_p = [p0, p_ref];
% state_v = [v0, v_ref];
% state_R = [R0, R_ref];
% state_w = [w0, w_ref];

%%
% load state_p.mat
% load state_v.mat
% load state_R.mat
% load state_w.mat
% load input_U.mat
% initial_input_U = input_U;
% initial_input_U(1,3)=initial_input_U(1,3)+100;
% initial_input_U(1,6)=initial_input_U(1,6)+100;
% initial_input_U(1,9)=initial_input_U(1,9)+100;
% initial_input_U(1,12)=initial_input_U(1,12)+100;

%%
objective = Objective_Function(state_p, state_v, state_R, state_w, initial_input_U, contact_sequence); obj_data(1) = objective;

tic;

prev_objective = 0;
initial_objective = objective;
minimum_objective = objective;
input_U = initial_input_U;
optimal_U = initial_input_U;
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

% if(max(max(abs(state_gap)))<1e-6)
%     step_size_alpha=1;
% end

    % Backward pass
    [del_u_feedforward, del_u_feedback_gain] = BackPass_DDP(input_U, contact_sequence);

    % Forward Pass
    [optimal_U, optimal_p, optimal_v, optimal_R, optimal_w] = ForwardPass_DDP(del_u_feedforward, del_u_feedback_gain, input_U);
    objective_forwardpass = Objective_Function([p0,optimal_p], [v0,optimal_v], [R0, optimal_R], [w0,optimal_w], optimal_U, contact_sequence);

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
    exit_condition1 = (max(max(abs(state_gap))) < epsilon_gap);% && (residual<1e-6);
    exit_condition2 = iter >= MAX_ITER;

    if ( exit_condition1 | exit_condition2 )
        break;
    end

end
fprintf("\n");
toc;

%%
state_error_psi(:,1) = Log_mat(R0);
for k=1:N
    state_error_psi(:,1+k) = Log_mat( R_ref(:,3*k-2:3*k).' * state_R(:,3*(k+1)-2:3*(k+1)) );
end

%% plots
figure; plot(0:iter,obj_data); grid;
title("Objective Function"); xlabel("iterations"); ylabel("Objective");

figure; plot(vecnorm(state_error_psi*180/pi));xlim([1,N+1]); grid; 
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
    rollout_R = R*Exp_mat(w*dt_mpc);
    rollout_w = w + w_dot*dt_mpc;

end


function objective = Objective_Function(state_p, state_v, state_R, state_w, input_U, contact_sequence)
    global weight_p weight_v weight_psi weight_w weight_u weight_cone N;
    global p_ref v_ref R_ref w_ref;
    global weight_state_bound x_lower x_upper;

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

    % nonlinear cone
    obj_cone = 0;
    for i = 1:N
        obj_cone = obj_cone + weight_cone.' * ( (input_U(:,i)-projection_cone_12x1(input_U(:,i),contact_sequence(:,i))).*(input_U(:,i)-projection_cone_12x1(input_U(:,i),contact_sequence(:,i))) );
    end

    predicted_x = [state_p(:,2:end);state_v(:,2:end);err_psi;state_w(:,2:end)];
    % state bounds
    obj_state_bound = 0;
    for i = 1:N
        bound_error = max(predicted_x(:,k)-x_upper, 0) + max(x_lower-predicted_x(:,k), 0);
        obj_state_bound = obj_state_bound + weight_state_bound.'*(bound_error.*bound_error);
    end

    objective = 0.5*(obj_p + obj_v + obj_psi + obj_w) + 0.5*obj_U + 0.5*obj_state_bound;% + 0.5*obj_cone;
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


function [y D] = projection_cone_3x1(force_3x1, contact)
    global mu_friction; global fz_max; global fz_min;

        y = zeros(3,1);
        D = zeros(3,3);

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
            D        = [1 0 0;
                        0 1 0;
                        0 0 0];
                
        elseif section2
            y(3,1)   = fz_max;
            y(1:2,1) = mu_friction*fz_max/radius * force_3x1(1:2);
            D        = mu_friction*fz_max/(radius^3)*...
                       [force_3x1(2)^2,  -force_3x1(1)*force_3x1(2),  0;
                       -force_3x1(1)*force_3x1(2),   force_3x1(1)^2,  0;
                                    0,                       0,           0];
            
        elseif section3
            y(3,1)   = (mu_friction*radius+force_3x1(3)) / (1+mu_friction^2);
            y(1:2,1) = mu_friction*y(3,1) * force_3x1(1:2)/radius;
            D        = [mu_friction^2/(1+mu_friction^2) + force_3x1(3)*force_3x1(2)^2/(radius^3),  -force_3x1(1)*force_3x1(2)*force_3x1(3)/(radius^3),  force_3x1(1)/radius;
                        -force_3x1(1)*force_3x1(2)*force_3x1(3)/(radius^3),  mu_friction^2/(1+mu_friction^2) + force_3x1(3)*force_3x1(1)^2/(radius^3),  force_3x1(2)/radius;
                        mu_friction/(1+mu_friction^2)*force_3x1(1)/radius,  mu_friction/(1+mu_friction^2)*force_3x1(2)/radius,  1/(1+mu_friction^2)];
                    
        elseif section4
            y(3,1)   = fz_min;
            y(1:2,1) = mu_friction*fz_min/radius * force_3x1(1:2);
            D        = mu_friction*fz_min/(radius^3)*...
                       [force_3x1(2)^2,  -force_3x1(1)*force_3x1(2),  0;
                       -force_3x1(1)*force_3x1(2),   force_3x1(1)^2,  0;
                            0,                       0,           0];
            
        elseif section5
            y(3,1)   = fz_min;
            y(1:2,1) = force_3x1(1:2);
            D        = [1 0 0;
                        0 1 0;
                        0 0 0];
            
        else
            y = force_3x1;
            D = eye(3);    
        end
    
    end
end



function [y D] = projection_cone_12x1(force12x1, contact_4leg)
    [y(1:3,1), D(1:3,:)]  = projection_cone_3x1(force12x1(1:3,1), contact_4leg(1));
    [y(4:6,1), D(4:6,:)]  = projection_cone_3x1(force12x1(4:6,1), contact_4leg(2));
    [y(7:9,1), D(7:9,:)]  = projection_cone_3x1(force12x1(7:9,1), contact_4leg(3));
    [y(10:12,1), D(10:12,:)] = projection_cone_3x1(force12x1(10:12,1), contact_4leg(4));

end


function grad_error = grad_state_bound(predicted_xk)
    global x_lower x_upper;
    
    grad_error = zeros(12,1);
    for i= 1:12
        if predicted_xk(i) > x_upper(i)
            grad_error(i) = 1;
        elseif predicted_xk(i) < x_lower(i)
            grad_error(i) = -1;
        else
            grad_error(i) = 0;
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
    global isProjected;
    global cone_mat dd_cone;
    global weight_state_bound x_lower x_upper;

    % initialization
    fx = zeros(12,12); 
    fu = zeros(12,12);
    I12x12 = eye(12);
    diff_cone_error = zeros(12,12);
    
    % weights
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

    % state_bounds
    state_x = [state_p;state_v;[init_err,state_err_psi];state_w];
    bound_error = max(state_x(:,N+1)-x_upper, 0) + max(x_lower-state_x(:,N+1), 0);
    grad_bound_error = grad_state_bound(state_x(:,N+1));
    jacobian_bound = weight_state_bound.*(grad_bound_error.*bound_error);
    hessian_bound = weight_state_bound.*(grad_bound_error.*grad_bound_error);

    Vxx = diag(terminal_weight) + diag(hessian_bound) + regularization_mu*I12x12;
    Vx  = terminal_weight .* state_err_x(:,N+1) + jacobian_bound + Vxx * state_gap(:,N+1);

    for backward_step = 1 : N
        k = (N+1) - backward_step;
        leg_position = init_leg_position - repmat(state_p(:,k)-p0,4,1);

        % nonlinear cone
        [projected_u D_cone] = projection_cone_12x1(input_U(:,k), contact_sequence(:,k));
        diff_cone_error(1:3,1:3) = I3 - D_cone(1:3,:);
        diff_cone_error(4:6,4:6) = I3 - D_cone(4:6,:);
        diff_cone_error(7:9,7:9) = I3 - D_cone(7:9,:);
        diff_cone_error(10:12,10:12) = I3 - D_cone(10:12,:);
        
        % state bounds
        bound_error = max(state_x(:,k)-x_upper, 0) + max(x_lower-state_x(:,k), 0);
        grad_bound_error = grad_state_bound(state_x(:,k));
        jacobian_bound = weight_state_bound.*(grad_bound_error.*bound_error);
        hessian_bound = weight_state_bound.*(grad_bound_error.*grad_bound_error);
        
        % force smooting
        % if backward_step==1
        %     err_u = zeros(12,1);
        % else
        %     err_u = input_U(:,k)-input_U(:,k+1);
        % end
        err_u = input_U(:,k);

        % model differentiation
        fx(1:3,1:3) = I3; fx(1:3,4:6) = I3*dt_mpc; 
        fx(4:6,4:6) = I3;
        fx(7:9,7:9)   = Exp_mat(state_w(:,k)*dt_mpc).';
        fx(7:9,10:12) = Right_Jacobian(state_w(:,k)*dt_mpc)*dt_mpc;
        fx(10:12,7:9)   = IB_INV * skew(state_R(:,3*k-2:3*k).'*skew_4legs(leg_position)*input_U(:,k)) * dt_mpc;
        fx(10:12,10:12) = I3 - IB_INV * ( skew(state_w(:,k))*IB - skew(IB*state_w(:,k)) ) * dt_mpc;

        fu(1:3,:) = 0.5*dt_mpc^2/m_mpc * I3x12;
        fu(4:6,:) = dt_mpc/m_mpc * I3x12;
        fu(10:12,:) = IB_INV * state_R(:,3*k-2:3*k).' * skew_4legs(leg_position) * dt_mpc;
        
        % Taylor approximation of Q-function
        Qx = weight_x.*state_err_x(:,k) + jacobian_bound + fx.'*Vx; 
        Qu = weight_u.*err_u + diff_cone_error.'*diag(weight_cone)*(input_U(:,k)-projected_u) + fu.'*Vx;
        Qxx = diag(weight_x) + diag(hessian_bound) + fx.'*Vxx*fx;
        Quu = diag(weight_u) + diff_cone_error.'*diag(weight_cone)*diff_cone_error + fu.'*Vxx*fu + regularization_mu*I12x12; Quu_INV = inv(Quu);
        Qux = fu.'*Vxx*fx;

        del_u_feedforward(:,k) = -Quu_INV*Qu;
        del_u_feedback_gain(:,12*k-11:12*k) = -Quu_INV*Qux;

        Vxx = Qxx + Qux.' * del_u_feedback_gain(:,12*k-11:12*k) + regularization_mu*I12x12; 
        Vx  = Qx  + Qux.' * del_u_feedforward(:,k) + Vxx*state_gap(:,k);  
        
        % Vx_data(:,k) = Vx; Vxx_data(:,12*k-11:12*k) = Vxx;
    end
    
end


function [optimal_U, optimal_p, optimal_v, optimal_R, optimal_w] = ForwardPass_DDP(del_u_feedforward, del_u_feedback_gain, input_U)
    global N;
    global p0 v0 R0 w0;
    global state_p state_v state_R state_w;
    global step_size_alpha;
    global state_gap;
    global delta1 delta2 Vx_data Vxx_data;

    % alpha step to the input direction
    del_x = zeros(12,1);

    optimal_U(:,1) = input_U(:,1) ...
                    + step_size_alpha * del_u_feedforward(:,1) ...
                    + del_u_feedback_gain(:,1:12)*del_x; 

    % alpha step to the state direction
    [rollout_p, rollout_v, rollout_R, rollout_w] = Euler_Integration([p0,v0,R0,w0], optimal_U(:,1));
    
    optimal_p(:,1)   = rollout_p + (step_size_alpha-1)*state_gap(1:3,1+1);
    optimal_v(:,1)   = rollout_v + (step_size_alpha-1)*state_gap(4:6,1+1);
    optimal_R(:,1:3) = rollout_R * Exp_mat((step_size_alpha-1)*state_gap(7:9,1+1));
    optimal_w(:,1)   = rollout_w + (step_size_alpha-1) * state_gap(10:12,1+1);

    for k = 2:N
        % alpha step to the input direction
        del_x = [optimal_p(:,k-1) - state_p(:,k);
                 optimal_v(:,k-1) - state_v(:,k);
                 Log_mat(state_R(:,3*(k)-2:3*(k)).'*optimal_R(:,3*(k-1)-2:3*(k-1))); 
                 optimal_w(:,k-1) - (state_R(:,3*(k)-2:3*(k)).'*optimal_R(:,3*(k-1)-2:3*(k-1))).'*state_w(:,k)];

        optimal_U(:,k) = input_U(:,k) ...
                + step_size_alpha * del_u_feedforward(:,k) ...
                + del_u_feedback_gain(:,12*k-11:12*k)*del_x;

% k
% (del_u_feedback_gain(:,12*k-11:12*k)*del_x).'
        
        % alpha step to the state direction
        [rollout_p, rollout_v, rollout_R, rollout_w] = Euler_Integration([optimal_p(:,k-1),optimal_v(:,k-1),optimal_R(:,3*(k-1)-2:3*(k-1)),optimal_w(:,k-1)], optimal_U(:,k));
        
        optimal_p(:,k)         = rollout_p + (step_size_alpha-1)*state_gap(1:3,1+k);
        optimal_v(:,k)         = rollout_v + (step_size_alpha-1)*state_gap(4:6,1+k);
        optimal_R(:,3*k-2:3*k) = rollout_R * Exp_mat((step_size_alpha-1)*state_gap(7:9,1+k));
        optimal_w(:,k)         = rollout_w + (step_size_alpha-1)*state_gap(10:12,1+k);

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

    % expected improvement
    % delta1 = delta1 + Qu_data(:,1).' * del_u_feedforward(:,1) ...
    %                 + state_gap(:,1).'*(Vx_data(:,1) - Vxx_data(:,1:12)*del_x);
    % delta2 = delta2 + del_u_feedforward(:,1).'*Quu_data(:,1:12)*del_u_feedforward(:,1) ...
    %                 + state_gap(:,1).'*(2*Vxx_data(:,1:12)*del_x - Vxx_data(:,1:12)*state_gap(:,1));


        % expected improvement 여기의 state_gap은 index +1 안하는게 맞다. 변수들 다 0~N에 맞는 인덱스로 되어 있음.
        % delta1 = delta1 + Qu_data(:,k).' * del_u_feedforward(:,k) ...
        %                 + state_gap(:,k).'*(Vx_data(:,k) - Vxx_data(:,12*k-11:12*k)*del_x);
        % delta2 = delta2 + del_u_feedforward(:,k).'*Quu_data(:,12*k-11:12*k)*del_u_feedforward(:,k) ...
        %                 + state_gap(:,k).'*(2*Vxx_data(:,12*k-11:12*k)*del_x - Vxx_data(:,12*k-11:12*k)*state_gap(:,k));
