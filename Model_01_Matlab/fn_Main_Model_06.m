
%Program: DLR-project 
%Topic: Motion and Parameter estimation of space objects using Laser-Vision
%data (Extended Kalman filter/Quaternion algebra)
%Author: Hrishik Mishra

%Function: fn_Main()
%Inputs: none
%Outputs:none
%Functionality: Initialize the system/Execute the Kalman filter/Display results
%12-state vector
%Implementing VB-EKF
%Author: Hrishik Mishra

function [ ] = fn_Main( )
    %Implementation of the Motion parameter estimation
    clc; close all;clear all;
    %Global variables
    global n;n = 0.0012; %Orbital velocity of Chaser/Target (nearly same)
    global parameter_gravitation;parameter_gravitation = 398.6005e12;
    global totalSimulationTime;totalSimulationTime = 50;     
    %VB-EKF specific details
    alpha = ones(6,1);
    beta = ones(6,1);
    dyn_param = 0.2;
    
    %%
     
        %Function: fn_Main()
        %Inputs: none
        %Outputs:none
        %Functionality: Read estimated pose file and display data
        %characteristics
        %Author: Hrishik Mishra
        %Measurements    
        %
    %1) Initialize the system with model parameters
    Init_Model_06;   
    
    %Integrate Dynamic equations to see result;    (Simulated data)
    %%
    dy_time = (0:t_delta:totalSimulationTime)';
    Time = dy_time;
    options = odeset('RelTol',1e-7,'AbsTol',1e-9*ones(13,1));
    X_a_0_sim = load('./Model_01_Matlab/X_a_0_Model_02');
    [T,X_a] = ode45(@fn_StateSpace,dy_time,X_a_0_sim.X_a_0_Model_02,options);
    save('Model_01_Matlab\X_a_true.mat','X_a');
    Signal_Quaternion = [];
    Mu = (fn_CrossTensor(ita,0)*X_a(:,1:4)')';
    Mu_noise = zeros(size(Mu));
    a = 50;
    b = 60;
%     Mu_noise(1:a,:) = (Mu(1:a,:) + 1e-4*randn(a,4));
%     Mu_noise(a+1:b,:) = (Mu(a+1:b,:) + 1e-4*randn(b-a,4));
%     Mu_noise(b+1:end,:) = (Mu(b+1:end,:) + 1e-4*randn(length(X_a(b+1:end,1)),4));   
    Mu_noise = quatnormalize(Mu+ normrnd(0,sqrt(3e-3),size(Mu)));
    r_c = zeros(3,length(dy_time));
    for iCount = 1:length(X_a)
        r_c(1:3,iCount) = X_a(iCount,8:10)' + rho_c + fn_CreateRotationMatrix(X_a(iCount,1:4)')*rho;
        Mu_noise(iCount,:) = Mu_noise(iCount,:)/norm(Mu_noise(iCount,:));
    end
    Signal_Quaternion(:,2:4) = Mu_noise(:,1:3);
    Signal_Quaternion(:,1) =Mu_noise(:,4); %+ %0.01*randn(length(X_a),1);
%     r_c(:,1:a) = r_c(:,1:a) + 3e-3*randn(a,3)';
%     r_c(:,a+1:b) = r_c(:,a+1:b) + 3e-3*randn(b-a,3)';
%     r_c(:,b+1:end) = r_c(:,b+1:end) + 3e-3*randn(length(X_a(b+1:end,1)),3)';
    r_c = r_c + normrnd(0,sqrt(3e-3),size(r_c));
    Signal_vector = timeseries(r_c(1:3,:),T,'Name','Signal');
    h_figure = figure('Name','Dynamic responses');
    subplot(2,2,1);
    plot(T,X_a(:,1),T,X_a(:,2),T,X_a(:,3),T,X_a(:,4),'LineWidth',2);
    legend('q_1', 'q_2', 'q_3','q_0');
    ylabel('Orientation');
    subplot(2,2,2)
    plot(T,X_a(:,5),T,X_a(:,6),T,X_a(:,7),'LineWidth',2);
    legend('\omega_x','\omega_y','\omega_z');
    ylabel('Angular rate');
    subplot(2,2,3);
    plot(T,X_a(:,8),T,X_a(:,9),T,X_a(:,10),'LineWidth',2);
    legend('r_x', 'r_y', 'r_z');
    ylabel('position(m)');
    subplot(2,2,4);
    plot(T,X_a(:,11),T,X_a(:,12),T,X_a(:,13),'LineWidth',2);
    legend('r_x_dot', 'r_y_dot', 'r_z_dot');
    ylabel('velocity(m/sec)');
    %%
    
    %Variables Initialization
    %%    
    X_pre = zeros(12, length(Time));%Time-series (12-state vector)        
    v_time = dy_time; %Time vector   
    X_post = zeros(12, length(Time) + 1);
    X_post(:,1) = X_0;
    
    X_a_Estimated = zeros(13,length(Time)+1);
    X_a_Estimated(:,1) = X_a_0;
    
    %Describe the state-error covariance initial matrix
    P_post = eye(12,12);    
%     P_post(1:3,1:3) = 0.5*eye(3,3);  
%     P_post(7:9,7:9) = 4*eye(3,3);
%     P_post(10:12,10:12) = eye(3,3);    
     
    residuals = zeros(length(v_time)-1,6);
    signal = zeros(7,1);
    test = zeros(6,length(v_time));
    options = odeset('RelTol',1e-7,'AbsTol',1e-9*ones(13,1));
    %%
    %%Kalman filter implementation in an iterative loop
    tic;
    for iCount = 2:length(v_time)
        %% Run the dynamic model to get predicted states        
        [~,X] = ode45(@fn_StateSpace,[v_time(iCount-1),v_time(iCount)],X_a_Estimated(:,iCount-1),options);                
        X_a_pre = X(end,:)';
        
        
        del_q = fn_CrossTensor(X_a_pre(1:4),0)*[-X_a_Estimated(1:3,iCount-1);X_a_Estimated(4,iCount-1)];
        a_g_pre = 2*del_q(1:3)/del_q(4);
        X_pre(:,iCount) = [a_g_pre;X_a_pre(5:13)];
        %% State propagation               
        Phi = fn_Create_Phi(X_post(:,iCount-1), n, t_delta,p);
        %% Create Process-noise covariance matrix
        Q_r = fn_Create_Q_r(Phi,X_pre(:,iCount), t_delta, sig_tau,sig_p,p);
        Q_t  = fn_Create_Q_t_k(t_delta, n,sig_f);
        Q = [Q_r, zeros(6,6);zeros(6,6),Q_t];
        
        
        %% Create State-error Covariance for Prediction-Stage
        P_pre = Phi*P_post*Phi' + Q;
        
        
        %% Fetch measurement
        signal(1:3,1) = Signal_vector.Data(:,:,iCount); 
        q_measured = Signal_Quaternion(iCount,:)';
        signal(4:7,1) = [q_measured(2:4); q_measured(1)];
        q_nominal = X_a_Estimated(1:4,iCount-1);
        %q_nominal = expm(0.5*fn_CrossTensor([X_a_Estimated(5:7,iCount-1);0],0)*t_delta)*X_a_Estimated(1:4,iCount-1);        
        zk =fn_Create_obs(signal,rho_c,q_nominal,ita);
        test(:,iCount-1) = zk;
        
        
        
        h = fn_Create_h(q_nominal,X_a_pre,rho);
        
        Hk = fn_Create_H(q_nominal,rho);
        
        Sk = fn_Create_S(q_nominal,ita,Cov_r,Cov_nu);
        
        K = fn_ComputeKalmanGain(P_pre,Hk,Sk);                    
        residuals(iCount-1,:) = zk - h;
        X_post(:,iCount) = X_pre(:,iCount) + K*(zk - h);
        P_post = (eye(12,12)-K*Hk)*P_pre;
        %Set predicted alpha/beta values
        
        alpha = dyn_param*alpha;
        beta = dyn_param*beta;
        betaprev = 0;
        alpha = alpha + 0.5;
        %VB iteration to converge to fixed values
        %for vbCounter = 1:10
        %% Update phase
            
            %Sk = fn_Create_S(q_nominal,ita,Cov_r,Cov_nu);
            %Sk = diag(beta./alpha);
            %K = fn_ComputeKalmanGain(P_pre,Hk,Sk);                    
            %residuals(iCount-1,:) = zk - h;
            %X_pre_estimate(:,iCount) = X_pre(:,end) + K*(zk - h);            
            %P_post = (eye(12,12)-K*Hk)*P_pre;
            %post_h = fn_Create_h(fn_CrossTensor([X_pre_estimate(1:3,iCount);1],0)*q_nominal,X_pre_estimate(:,iCount), rho);
            %post_residual = zk - post_h;
            %matA = Hk*P_post*Hk';
            %beta = beta + 0.5*post_residual.^2 + 0.5*diag(matA);
            %if norm(beta - betaprev) < 1e-6
%                fprintf('counter:%d\n',vbCounter);
%                break;                
%            end
%            betaprev = beta;
%        end
        
       %Convert to quaternion
       a_g = X_post(1:3,iCount);
       %rho_q = fn_CrossTensor([a_g;2],0)*q_nominal;
       del_q = [a_g;2]/(sqrt(4+dot(a_g,a_g)));
       %q = rho_q/norm(rho_q); 
       q = fn_CrossTensor(del_q,0)*q_nominal;
       X_a_Estimated(1:4,iCount) = q;
       X_a_Estimated(5:13,iCount) = X_post(4:12,iCount);
       %Reset quat vector
       X_post(1:3,iCount) = 0;
        
    end
    X_a_Estimated = X_a_Estimated(:,2:end);
    
     %% Generate measurements to compare with actual measurements
    Mu_est = zeros(4,length(X_a_Estimated));
    r_c_est = zeros(3,length(X_a_Estimated));
    for iCount = 1:length(X_a_Estimated)
       Mu_est(:,iCount) = fn_CrossTensor(ita,0)*X_a_Estimated(1:4,iCount);
       r_c_est(1:3,iCount) = X_a_Estimated(8:10,iCount) + rho_c + fn_CreateRotationMatrix(X_a_Estimated(1:4,iCount))*rho;
    end
    
    %% Plot the details
    figure('Name','Post processing details');
    subplot(4,2,1);
    stairs(dy_time,reshape(Signal_vector.Data(1,:,:),length(Signal_vector.Data(1,:,:)),1),'-');
    hold all;
    plot(dy_time,r_c_est(1,:),'LineWidth',2);
    legend('Measured','Estimated');
    ylabel('x-coordinate');
    subplot(4,2,2);
    stairs(dy_time,reshape(Signal_vector.Data(2,:,:),length(Signal_vector.Data(2,:,:)),1),'-');
    hold all;
    plot(dy_time,r_c_est(2,:),'LineWidth',2);
    legend('Measured','Estimated');
    ylabel('y-coordinate');
    subplot(4,2,3);
    stairs(dy_time,reshape(Signal_vector.Data(3,:,:),length(Signal_vector.Data(3,:,:)),1),'-');
    hold all;
    plot(dy_time,r_c_est(3,:),'LineWidth',2);
    legend('Measured','Estimated');
    ylabel('z-coordinate');
    
    subplot(4,2,5);
    stairs(dy_time,Signal_Quaternion(:,1));
    hold all;
    plot(dy_time,Mu_est(4,:),'LineWidth',2);
    legend('Measured','Estimated');
    ylabel('q_0');
    subplot(4,2,6);
    stairs(dy_time,Signal_Quaternion(:,2));
    hold all;
    plot(dy_time,Mu_est(1,:),'LineWidth',2);
    legend('Measured','Estimated');
    ylabel('q_1');
    subplot(4,2,7);
    stairs(dy_time,Signal_Quaternion(:,3));
    hold all;
    plot(dy_time,Mu_est(2,:),'LineWidth',2);
    legend('Measured','Estimated');
    ylabel('q_2');
    subplot(4,2,8);
    stairs(dy_time,Signal_Quaternion(:,4));
    hold all;
    plot(dy_time,Mu_est(3,:),'LineWidth',2);
    legend('Measured','Estimated');
    ylabel('q_3');
    
    figure;
    subplot(2,1,1);
    plot(dy_time,X_a_Estimated(5:7,:)');hold all;
    plot(dy_time,X_a(:,5:7),'LineStyle','-.');
    legend('x_est','y_est','z_est','x_true','y_true','z_true');
    ylabel('\omega');
    
%     subplot(3,1,2);
%     plot(dy_time,X_a_Estimated(8:10,:)');hold all;  
%     plot(dy_time,X_a(:,8:10),'LineStyle','-.')
%     legend('x_est','y_est','z_est','x_true','y_true','z_true');
%     ylabel('p');
    subplot(2,1,2);
    plot(dy_time,X_a_Estimated(11:13,:)');hold all;
    plot(dy_time,X_a(:,11:13),'LineStyle','-.');
    legend('v_x_est','v_y_est','v_z_est','v_x_true','v_y_true','v_z_true');
    ylabel('velocity');    
    
    
    figure;
    subplot(3,2,1);
    plot(dy_time,X_a_Estimated(1,:)');hold all;
    plot(dy_time,X_a(:,1),'LineStyle','-.');
    legend('q_1 est','q_1 true');
    ylabel('quaternion');
    subplot(3,2,2);    
    plot(dy_time,X_a_Estimated(2,:)');hold all;
    plot(dy_time,X_a(:,2),'LineStyle','-.');
    legend('q_2 est','q_2 true');
    ylabel('quaternion');
    
    subplot(3,2,3);
    plot(dy_time,X_a_Estimated(3,:)');hold all;
    plot(dy_time,X_a(:,3),'LineStyle','-.');
    legend('q_3 est','q_3 true');
    ylabel('quaternion');
    
    subplot(3,2,4);
    
    plot(dy_time,X_a_Estimated(4,:)');hold all;
    plot(dy_time,X_a(:,4),'LineStyle','-.');
    legend('q_0 est','q_0 true');
    ylabel('quaternion');
    
    subplot(3,2,5);
    plot(dy_time,X_a_Estimated(11:13,:));hold all;
    plot(dy_time,X_a(:,11:13),'LineStyle','-.');
    legend('r_x est','r_y est', 'r_z est', 'r_x true','r_y true', 'r_z true');
    
    subplot(3,2,6);
    plot(dy_time(2:end),residuals(:,4:6)');hold all;
    %legend('r_x est','r_y est', 'r_z est', 'r_x true','r_y true', 'r_z true');
    
    save('Model_01_Matlab\X_a_Estimated_2.mat','X_a_Estimated');
    
end

%Function: fn_StateSpace()
%Inputs: X_a (23-state vector)
%Outputs: [dy] (Time-derivative of the state vector)
%Functionality: Implements the State-space model of the dynamic system
%Author: Hrishik Mishra
function[dy] = fn_StateSpace(~,X_a)
    global n;
    Init_Model_06;
    v_n = [0;0; n];
    q = X_a(1:4);
    omega = X_a(5:7);
    p = [p_x;p_y;p_z];
    r = X_a(8:10);
    r_dot = X_a(11:13);
    rho_t = rho;
    ita_t = ita;
    
    J_k = [1 0 0;0 (1-p(2))/(1+p(1)) 0;0 0 (1+p(3))/(1-p(1))];
    q_omega_rel = [omega;0];%Check this for the Attitude Control of the Chaser
    psi = [p(1)*omega(2)*omega(3); p(2)*omega(1)*omega(3);p(3)*omega(1)*omega(2)];
    q_dot = 0.5*fn_CrossTensor(q_omega_rel,0)*q;
    tau = [tau_1;tau_2;tau_3]*rand;
    e_force = [e_force_x;e_force_y;e_force_z]*rand;
    omega_dot = psi + J_k*tau.*rand(3,1); %for static conditions
    rho_t_dot = zeros(3,1);
    r_ddot = fn_Compute_r_ddot(n,r,e_force,r_dot);
    %r_dot = zeros(3,1);
    %r_ddot = zeros(3,1);
    dy = [q_dot;omega_dot;r_dot;r_ddot];
end


%Function: fn_Compute_r_ddot()
%Inputs: n - angular velocity of spacecrafts,
%        r - position vector between Chaser CoM and Target CoM
%        e_force - 3-D force vector for perturbation
%        r_dot - velocity-vector between Chaser CoM and Target CoM
%Outputs: r_ddot (Time-derivative of the velocity vector)
%Functionality: Computes time-derivative of velocity-vector between Chaser
%CoM and Target CoM
%Author: Hrishik Mishra
function r_ddot = fn_Compute_r_ddot(n,r,e_force,r_dot)
%
    global parameter_gravitation;
    re = [(parameter_gravitation/n^2)^(1/3);0;0];
    v_n = [0;0;n];
    Tensor_n = fn_VectorToSkewSymmetricTensor(v_n);
    Term_1 = -2*Tensor_n*r_dot;
    Term_2 = Tensor_n*(Tensor_n*r);
    Term_3 = parameter_gravitation*(re + r)/(norm(re+r))^3;
    Term_4 = (n^2)*re;
    
    r_ddot = Term_1 - Term_2 - Term_3 + Term_4 + e_force;
end

%Function: fn_CreateRotationMatrix()
%Inputs: q: quaternion between Chaser CoM and Target CoM
%Outputs: R (Rotation Matrix [3,3])
%Functionality: Generates R
%Author: Hrishik Mishra
function [R] = fn_CreateRotationMatrix(q)
    q_0 = q(4);
    q_v = q(1:3);
    R = (2*q_0^2-1)*eye(3,3) + 2*q_0*fn_VectorToSkewSymmetricTensor(q_v) + 2*q_v*q_v';
end

%Create state propagation matrix of the whole system

%Function: fn_Create_Phi()
%Inputs: X_a - 23-state vector,
%        n - angular velocity of spacecrafts,
%        t_delta - sampling time
%Outputs: [Phi] (State-propagation Matrix)
%Functionality: Generates the State Propagation matrix
%Author: Hrishik Mishra
function[Phi] = fn_Create_Phi(X,n,t_delta,p)
    Phi_t = fn_Create_Phi_t(n,t_delta);
    Phi_r = fn_Create_Phi_r(X,p);
    Phi = [Phi_r, zeros(6,6);zeros(6,6),Phi_t];

end
%Create state propagation matrix of the whole system

%Function: fn_Create_Phi_r()
%Inputs: X_a - 23-state vector,
%Outputs: [Phi_r] (State-propagation Matrix for Rotation component)
%Functionality: Generates the State Propagation matrix for Rotation
%component
%Author: Hrishik Mishra
function [Phi_r] = fn_Create_Phi_r(X,p)
    global t_delta;    
    omega = X(4:6);
    M = fn_Create_M(p,omega);   
    
    A = [-fn_VectorToSkewSymmetricTensor(omega),eye(3,3);zeros(3,3),M];
    Phi_r = expm(A*t_delta);
%     phi_r11 = fn_Create_phi_r11(omega,t_delta);
%     phi_r12 = fn_Create_phi_r12(omega,M,t_delta);
%     phi_r22 = fn_Create_phi_r22(omega,t_delta,M);
%     phi_r13 = fn_Create_phi_r13(omega,t_delta,M,N);
%     phi_r23 = fn_Create_phi_r23(omega,t_delta,M,N);
%     Phi_r = [phi_r11, phi_r12, phi_r13;zeros(3,3),phi_r22,phi_r23;zeros(3,3),zeros(3,3),eye(3,3)];
%     
end
%Function: fn_Create_M()
%Inputs: p - Inertial ratios 
%        omega - angular velocity between Chaser Grasping frame and Chaser
%        CoM
%Outputs: M (Linearized model matrix for inertial parameters)
%Functionality: Generates M
%Author: Hrishik Mishra
function M = fn_Create_M(p,omega)
%
    M = [0, p(1)*omega(3), p(1)*omega(2); p(2)*omega(3), 0, p(2)*omega(1); p(3)*omega(2),p(3)*omega(1),0];
end


%%
%Function: fn_Create_Phi_t()
%Inputs: n - angular velocity of Chaser/Target (almost same)
%        t_delta - sampling time
%Outputs: Phi_t (State Transition Matrix for Translation component)
%Functionality: Generates Phi_t
%Author: Hrishik Mishra
function Phi_t = fn_Create_Phi_t(n,t_delta)
%
    Phi_t12 = [t_delta, n*t_delta^2,0;0, t_delta, 0;0,0,t_delta];
    Phi_t22 = eye(3,3) - 2*fn_VectorToSkewSymmetricTensor([0;0;n])*t_delta;
    
    Phi_t = [eye(3,3),Phi_t12;zeros(3,3),Phi_t22];
    
end

%Function: fn_Create_Q_t_k()
%Inputs: t_delta: sampling time
%        sigma_f: variance of force perturbations
%Outputs: Q_t_k (Translational component Process noise)
%Functionality: Generates Q_t_k
%Author: Hrishik Mishra
function Q_t_k  = fn_Create_Q_t_k(t_delta, n,sigma_f)
%

    Q_t11_k = [(t_delta^3)/3 + (2/5)*n^2*t_delta^5, (1/4)*(n - n^2)*t_delta^4, 0;(1/4)*(n - n^2)*t_delta^4, (t_delta^3)/3 - (4/15)*n^2*t_delta^5 0;0, 0, (t_delta^3)/3 - (1/15)*n^2*t_delta^5];
    Q_t12_k = [(t_delta^2)/2 + (1/3)*(n^2)*(t_delta^4), -(1/3)*n*(t_delta^3),0;(2/3)*(n - n^2)*t_delta^3, (1/2)*(t_delta^2) - (2/3)*(n^2)*(t_delta^4) 0;0 0 (1/2)*(t_delta^2) - (1/6)*(n^2)*(t_delta^4)];
    Q_t22_k = [t_delta + (n^2)*(t_delta^3), 0, 0;0, t_delta, 0;0, 0, t_delta - (1/3)*(n^2)*(t_delta^3)];

    Q_t_k = [Q_t11_k, Q_t12_k;Q_t12_k',Q_t22_k];

    Q_t_k = sigma_f*Q_t_k;
end


%Function: fn_Create_Q_r()
%Inputs: Phi - State Propagation Matrix for Rotational component
%        X - 23-state vector
%        t_delta - sampling time
%        tau - Torque disturbance
%        sig_p - Uncertainty in Inertial ratios (p)
%Outputs: Q_r (Rotation component Process noise)
%Functionality: Generates Q_r
%Author: Hrishik Mishra
function Q_r = fn_Create_Q_r(Phi,X, t_delta, tau,sig_p,p)
    J = fn_Create_J(p);
    Q_r11 = 0.25*fn_Create_Q_r11(Phi(1:3,4:6), J, t_delta, tau);
    Q_r12 = 0.5*fn_Create_Q_r12(Phi(1:3,4:6),Phi(4:6,4:6),J,t_delta,tau);
    Q_r22 = fn_Create_Q_r22(Phi(4:6,4:6),J,t_delta,tau);
    Q_r = [Q_r11, Q_r12;Q_r12' Q_r22];
    %check if it is actually phi_r12 or phi_r11
    function Q_r11 = fn_Create_Q_r11(phi_r12, J, t_delta, tau)
        Q_r11 = tau*phi_r12*J^2*phi_r12'*t_delta; 
    end

    function Q_r12 = fn_Create_Q_r12(phi_r12,phi_r22,J, t_delta, tau)
        Q_r12 = tau*phi_r12*J^2*phi_r22*t_delta;
    end
    function Q_r22 = fn_Create_Q_r22(phi_r22, J, t_delta, tau)
        Q_r22 = tau*phi_r22*J^2*phi_r22*t_delta;
    end
end


%Function: fn_Create_obs()
%Inputs: signal: 7x1 measurement vector
%        rho_c: position vector of camera from Chaser CoM
%        q_k: quaternion between Chaser CoM and Target CoM
%        ita_k: quaternion between Target CoM and Grasping frame
%Outputs: zk (observation in terms of state-vector)
%Functionality: Generates zk
%Author: Hrishik Mishra
function zk =fn_Create_obs(signal,rho_c,q_k,ita_k)
%
    zk = zeros(6,1);
    mu = signal(4:7);
    rc = signal(1:3);
    
    ita_star = [-ita_k(1:3);ita_k(4)];
    q_star = [-q_k(1:3);q_k(4)];
    zk(1:3,1) = rc - rho_c;
    v = fn_CrossTensor(ita_star,0)*fn_CrossTensor(mu,0)*q_star;
    zk(4:6,1) = 2*v(1:3)/v(4);
end



%Function: fn_Create_h()
%Inputs: q: quaternion between Chaser CoM and Target CoM
%        X: 21-state vector 
%Outputs: h (non-linear measurement function)
%Functionality: Generates h
%Author: Hrishik Mishra
function h = fn_Create_h(q_nominal,X_a,rho)
%
    r = X_a(8:10,1);
    rho_t = rho;
    q = X_a(1:4);            
    del_q = fn_CrossTensor(q,0)*[-q_nominal(1:3);q_nominal(4)];
    a_g = 2*del_q(1:3)/del_q(4);    
    q_0 = q(4);
    q_v = q(1:3);
    Q_v = fn_VectorToSkewSymmetricTensor(q_v);
    R_q = (2*q_0^2 - 1)*eye(3,3) + 2*q_0*Q_v + 2*(q_v)*(q_v)';
    h1 = r + R_q*rho_t;
    temp = a_g;
    h2 = temp(1:3);
    h = [h1;h2];    
end


%Function: fn_Create_S()
%Inputs: q_nominal: Quaternion between Chaser CoM and Target CoM
%        ita_nominal: Quaternion between Target CoM and Target grasping
%        frame
%        cov_r: position covariance due to camera
%        cov_nu: orientation covariance due to camera
%Outputs: S (Measurement noise covariance Matrix)
%Functionality: Generates S
%Author: Hrishik Mishra
function S   = fn_Create_S(q_nominal,ita_nominal,cov_r,cov_nu)
%
    q_star = [-q_nominal(1:3);q_nominal(4)];
    ita_star = [-ita_nominal(1:3);ita_nominal(4)];
    Ita_Tensor = fn_CrossTensor(ita_star,0);
    Q_Tensor = fn_CrossTensor(q_star,1);
    T = [eye(3,3),zeros(3,1)]*Ita_Tensor*Q_Tensor;
    S = [cov_r,zeros(3,3);zeros(3,3), T*cov_nu*T'];
end
%Function: fn_ComputeKalmanGain()
%Inputs: Pk - Predicted state-error Covariance Matrix
%        Hk - Linearized measurement function Matrix
%        Sk - Measurement covariance Matrix
%Outputs: K (Kalman gain)
%Functionality: Generates K
%Author: Hrishik Mishra
function K = fn_ComputeKalmanGain(Pk,Hk,Sk)
%
    K = Pk*Hk'/(Hk*Pk*Hk' + Sk);
end



%Function: fn_Create_H()
%Inputs: q_nominal: Quaternion between Chaser CoM and Target CoM
%        X: 21-state vector
%Outputs: H (Linearized measurement function matrix)
%Functionality: Generates H
%Author: Hrishik Mishra
function H = fn_Create_H(q_nominal, rho)    
    Matrix_1 = fn_Create_del_h1(q_nominal,rho);       
    H = [Matrix_1;eye(3,3), zeros(3,9)];
    function Matrix_1 = fn_Create_del_h1(q,rho_t_k)
        q_0 = q(4);
        q_v = q(1:3);
        Q_v = fn_VectorToSkewSymmetricTensor(q_v);
        R = (2*q_0^2 - 1)*eye(3,3) + 2*q_0*Q_v + 2*(q_v)*(q_v)';

        first = -2*R*fn_VectorToSkewSymmetricTensor(rho_t_k);
        Matrix_1 = [first, zeros(3,3),eye(3,3),zeros(3,3)];
    end    
end

%Function: fn_Create_J()
%Inputs: p: Inertial ratios
%Outputs: J (Inertial Ratio Matrix for omega_dot)
%Functionality: Generates J
%Author: Hrishik Mishra
function J = fn_Create_J(p)
    J = [1 0 0;0 (1-p(2))/(1+p(1)) 0;0 0 (1+p(3))/(1-p(1))];
end
