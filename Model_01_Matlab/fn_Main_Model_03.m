
%Program: DLR-project 
%Topic: Motion and Parameter estimation of space objects using Laser-Vision
%data (Extended Kalman filter/Quaternion algebra)
%Author: Hrishik Mishra

%Function: fn_Main()
%Inputs: none
%Outputs:none
%Functionality: Initialize the system/Execute the Kalman filter/Display results
%Author: Hrishik Mishra

function [ ] = fn_Main_Model_03( )
    %Implementation of the Motion parameter estimation
    clc; close all;clear all;
    %Global variables
    global n;n = 0.0012; %Orbital velocity of Chaser/Target (nearly same)
    global parameter_gravitation;parameter_gravitation = 398.6005e12;
    global totalSimulationTime;totalSimulationTime = 180;        
    
    %%
     
        %Function: fn_Main()
        %Inputs: none
        %Outputs:none
        %Functionality: Read estimated pose file and display data
        %characteristics
        %Author: Hrishik Mishra
        %Measurements    
        (read_estimates_pose); %Read the Pose data from the file provided
        %Signal_Transf = Transformation (DCM) between Camera-frame and
        %Grasping-frame
        Signal_Transf = ts_coll.Signal_Transf_2;
        %Signal_vector = Position vector between Camera-frame and Grasping-
        %frame
        Signal_vector = ts_coll.Signal_vector_2;
        Time = Signal_Transf.Time; Time = Time + 0.1;

        h_figure_1 = figure('Name','Measurement details');
        subplot(3,1,1);
        h_stem = stem(Time,reshape(Signal_vector.Data(1,1,:),length(Time),1),'Marker','o','LineStyle','-.','MarkerFaceColor','red');
        set(h_stem,'BaseValue',mean(reshape(Signal_vector.Data(1,1,:),length(Time),1)));
        ylabel('x-measurement frame');
        subplot(3,1,2);
        h_stem = stem(Time,reshape(Signal_vector.Data(2,1,:),length(Time),1),'Marker','o','LineStyle','-.','MarkerFaceColor','green');
        set(h_stem,'BaseValue',mean(reshape(Signal_vector.Data(2,1,:),length(Time),1)));
        ylabel('y-measurement frame');
        subplot(3,1,3);
        h_stem = stem(Time,reshape(Signal_vector.Data(3,1,:),length(Time),1),'Marker','o','LineStyle','-.','MarkerFaceColor','green');
        set(h_stem,'BaseValue',mean(reshape(Signal_vector.Data(3,1,:),length(Time),1)));
        ylabel('z-measurement frame');

        %Convert measurement orientation to quaternions
        Signal_Quaternion = dcm2quat(Signal_Transf.Data);
        h_figure = figure('Name','Orientation details');
        subplot(2,1,1);
        h1_stem = stem(Time,quatnorm(dcm2quat(Signal_Transf.Data)),'LineStyle','-.','Marker','o','MarkerFaceColor','red');
        set(h1_stem, 'BaseValue',mean(quatnorm(dcm2quat(Signal_Transf.Data))));
    
    
    %1) Initialize the system with model parameters
    Init_Model_01_Matlab;   
    
    %Integrate Dynamic equations to see result;    (Simulated data)
    %%
    dy_time = (0:t_delta:totalSimulationTime)';
    Time = dy_time;
    options = odeset('RelTol',1e-7,'AbsTol',1e-9*ones(23,1));
    X_a_0_sim = load('./Model_01_Matlab/X_a_0_est');
    [T,X_a] = ode45(@fn_StateSpace,dy_time,X_a_0_sim.X_a_0_sim,options);
    Signal_Quaternion = [];
    Mu = (fn_CrossTensor(X_a(1,20:23)',0)*X_a(:,1:4)')';
    Mu_noise = quatnormalize(Mu + 5e-2*randn(length(X_a),4));
    Signal_Quaternion(:,2:4) = Mu_noise(:,1:3);
    Signal_Quaternion(:,1) = X_a(:,4); %+ %0.01*randn(length(X_a),1);
    r_c = zeros(3,length(dy_time));
    for iCount = 1:length(X_a)
        r_c(1:3,iCount) = X_a(iCount,11:13)' + rho_c + fn_CreateRotationMatrix(X_a(iCount,1:4)')*X_a(iCount,17:19)';
    end
    Signal_vector = timeseries(r_c(1:3,:) + 3e-2*randn(length(r_c),3)',T,'Name','Signal');
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
    plot(T,X_a(:,11),T,X_a(:,12),T,X_a(:,13),'LineWidth',2);
    legend('r_x', 'r_y', 'r_z');
    ylabel('position(m)');
    subplot(2,2,4);
    plot(T,X_a(:,14),T,X_a(:,15),T,X_a(:,16),'LineWidth',2);
    legend('r_x_dot', 'r_y_dot', 'r_z_dot');
    ylabel('velocity(m/sec)');
    %%
    
    %Variables Initialization
    %%
    X_a_init = X_a_0; %23-state vector to store the final output of Integration of dynamic equations (Predicted Result)
    X_a_pre = zeros(23, length(Time));%Time-series of X_a_init (23-state vector)
    X_a_itr = X_a_0;%23-state vector to store the output after Update-Stage (Estimated result)    
    X_a_Estimated = zeros(23, length(Time));%Time-series of X_a_itr
    X_a_Estimated(:,1) = X_a_itr;%Initialize the vector
    X_pre = zeros(21,1);%21-state vector to store the output after Prediction-Stage
    X_pre_estimate = zeros(21,length(Time));%Time-series of 21-state vector after Update-stage    
    X_a_pre(:,1) = X_a_init;
    v_time = dy_time; %Time vector   
    X_pre(4:18,1) = X_a_0(5:19,1);
    
    
    %Describe the state-error covariance initial matrix
    P_post = eye(21,21);    
%     P_post(1:3,1:3) = 0.5*eye(3,3);
%     P_post(7:9,7:9) = 0.005*eye(3,3);   
%     P_post(10:12,10:12) = 4*eye(3,3);
%     P_post(13:15,13:15) = eye(3,3);
%     P_post(16:18,16:18) = 0.04*eye(3,3);
%     P_post(19:21,19:21) = 0.005*eye(3,3);
     
    residuals = zeros(length(v_time)-1,6);
    signal = zeros(7,1);
    %%
    %%Kalman filter implementation in an iterative loop
    tic;
    for iCount = 2:length(v_time)
        %q_nominal = nominal-quaternion (Chaser CoM-Target CoM) from previous states
        q_nominal = expm(0.5*fn_CrossTensor([X_a_itr(5:7,1);0],0)*t_delta)*X_a_itr(1:4,1);
        %ita_nominal = nominal-quaternion (Target CoM-Target Frame) from
        %previous states
        ita_nominal = X_a_itr(20:23,1);
        
        %% Run the dynamic model to get predicted states        
        [~,X] = ode45(@fn_StateSpace,[v_time(iCount-1),v_time(iCount)],X_a_itr,options);                
        X_a_init = X(end,:)';
        X_a_pre(:,iCount) = X_a_init;
        
                
        %% Change to intermediate state
        del_q = fn_CrossTensor(X_a_init(1:4,end),0)*[-q_nominal(1:3);q_nominal(4)];
        ita_star = [-ita_nominal(1:3);ita_nominal(4)];
        del_ita = fn_CrossTensor(ita_star,0)*X_a_init(20:23,end);
        %Set up X_pre for the Update-stage
        X_pre(1:3,end) = del_q(1:3); 
        X_pre(4:18,end) = X_a_init(5:19,end);
        X_pre(19:21,end) = del_ita(1:3);
        
        
        %% State propagation               
        Phi = fn_Create_Phi(X_a_itr, n, t_delta);
        
        
        %% Create Process-noise covariance matrix
        Q_r = fn_Create_Q_r(Phi,X_a_itr, t_delta, sig_tau,sig_p);
        Q_t  = fn_Create_Q_t_k(t_delta, n,sig_f);
        Q_theta = sig_theta*eye(6,6);
        Q = [Q_r, zeros(9,6), zeros(9,6);zeros(6,9), Q_t, zeros(6,6);zeros(6,9),zeros(6,6),Q_theta];

        
        %% Create State-error Covariance for Prediction-Stage
        P_pre = Phi*P_post*Phi' + Q;
        
        
        %% Fetch measurement
        signal(1:3,1) = Signal_vector.Data(:,:,iCount); 
        q_measured = Signal_Quaternion(iCount,:)';
        signal(4:7,1) = [q_measured(2:4); q_measured(1)];
        zk =fn_Create_obs(signal,rho_c,q_nominal,ita_nominal);
        
        
        %% Update phase
        Hk = fn_Create_H(q_nominal,X_pre(:,end));
        Sk = fn_Create_S(q_nominal,ita_nominal,Cov_r,Cov_nu);
        K = fn_ComputeKalmanGain(P_pre,Hk,Sk);
        h = fn_Create_h(fn_CrossTensor(del_q,0)*q_nominal,X_pre);        
        residuals(iCount-1,:) = zk - h;
        X_pre_estimate(:,iCount) = X_pre(:,end) + K*(zk - h);
        
        
        %% Error/Warning checks in computed Quaternions
        del_q_v = X_pre_estimate(1:3,iCount);
        if ( norm(del_q_v) > 1)
           display('Warning: Normalizing del_q_v');
           del_q_v = del_q_v/norm(del_q_v);
           del_q_0 = 0;
        else
            del_q_0 = sqrt(1 - norm(del_q_v)^2);
        end        
        del_q = [del_q_v;del_q_0];
        
        if (imag(del_q(1)) ~= 0)
            display('Error due to del_q being imaginary');
            break;
        end
        if ( abs(del_q_0) < 0.8 )
            display('Warning, del_q_0 is too low');
            fprintf('del_q_0:%f\n',del_q_0);
        end        
       
         
        del_ita_v = X_pre_estimate(19:21,iCount);
        if ( norm(del_ita_v) > 1)
           display('Warning: Normalizing del_ita_v');
           del_ita_v = del_ita_v/norm(del_ita_v); 
           del_ita_0 = 0;
        else
            del_ita_0 = sqrt(1 - norm(del_ita_v)^2);
        end
        
        del_ita = [del_ita_v;del_ita_0];
        
        if (imag(del_ita(1))~= 0)
            display('Error due to del_ita being imaginary');
            break;
        end
        if ( abs(del_ita_0) < 0.7 )
            display('Warning, del_ita_0 is too low');
            fprintf('del_ita_0:%f\n',del_ita_0);
        end
        
         
        %% Convert from intermediate to Estimated-states after Update-stage 
        q_est = fn_CrossTensor(del_q,0)*q_nominal;
        ita_est = fn_CrossTensor(del_ita,1)*ita_nominal;        
        X_a_itr = [q_est;X_pre_estimate(4:18,iCount);ita_est];
        X_a_Estimated(:,iCount) = X_a_itr;
        P_post = (eye(21,21)-K*Hk)*P_pre;
  
        
    end
    toc;
    %% Generate measurements to compare with actual measurements
    Mu_est = zeros(4,length(X_a_Estimated));
    r_c_est = zeros(3,length(X_a_Estimated));
    for iCount = 1:length(X_a_Estimated)
       Mu_est(:,iCount) = fn_CrossTensor(X_a_Estimated(20:23,iCount),0)*X_a_Estimated(1:4,iCount);
       r_c_est(1:3,iCount) = X_a_Estimated(11:13,iCount) + rho_c + fn_CreateRotationMatrix(X_a_Estimated(1:4,iCount))*X_a_Estimated(17:19,iCount);
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
    subplot(3,1,1);
    plot(dy_time,X_a_Estimated(5:7,:)');hold all;
    plot(dy_time,X_a(:,5:7),'LineStyle','-.');
    legend('x_est','y_est','z_est','x_true','y_true','z_true');
    ylabel('\omega');
    
    subplot(3,1,2);
    plot(dy_time,X_a_Estimated(8:10,:)');hold all;  
    plot(dy_time,X_a(:,8:10),'LineStyle','-.')
    legend('x_est','y_est','z_est','x_true','y_true','z_true');
    ylabel('p');
    subplot(3,1,3);
    plot(dy_time,X_a_Estimated(14:16,:)');hold all;
    plot(dy_time,X_a(:,14:16),'LineStyle','-.');
    legend('v_x_est','v_y_est','v_z_est','v_x_true','v_y_true','v_z_true');
    ylabel('velocity');
    
    figure;
    subplot(2,1,1);
    plot(dy_time,X_a_Estimated(17:19,:)');hold all;
    plot(dy_time,X_a(:,17:19),'LineStyle','-.');
    legend('\rho_tx_est','\rho_ty_est','\rho_tz_est','\rho_tx_true','\rho_ty_true','\rho_tz_true');
    ylabel('\rho_t');
    subplot(2,1,2);
    plot(dy_time,X_a_Estimated(20:23,:)');hold all;
    plot(dy_time,X_a(:,20:23),'LineStyle','-.');
    legend('\eta_1_est','\eta_2_est','\eta_3_est','\eta_0_est','\eta_1_true','\eta_2_true','\eta_3_true','\eta_0_true');
    ylabel('\eta');
    
    figure;
    subplot(3,1,1);
    plot(dy_time,X_a_Estimated(1:4,:)');hold all;
    plot(dy_time,X_a(:,1:4),'LineStyle','-.');
    legend('q_1 est','q_2 est','q_3 est','q_0 est','q_1 true','q_2 true','q_3 true','q_0 true');
    ylabel('quaternion');
    subplot(3,1,2);
    plot(dy_time,X_a_Estimated(11:13,:));hold all;
    plot(dy_time,X_a(:,11:13),'LineStyle','-.');
    legend('r_x est','r_y est', 'r_z est', 'r_x true','r_y true', 'r_z true');
    
    
end

%Function: fn_StateSpace()
%Inputs: X_a (23-state vector)
%Outputs: [dy] (Time-derivative of the state vector)
%Functionality: Implements the State-space model of the dynamic system
%Author: Hrishik Mishra
function[dy] = fn_StateSpace(~,X_a)
    global n;
    Init_Model_01_Matlab;
    v_n = [0;0; n];
    q = X_a(1:4);
    omega = X_a(5:7);
    p = X_a(8:10);
    r = X_a(11:13);
    r_dot = X_a(14:16);
    rho_t = X_a(17:19);
    ita_t = X_a(20:23);
    
    J_k = [1 0 0;0 (1-p(2))/(1+p(1)) 0;0 0 (1+p(3))/(1-p(1))];
    q_omega_rel = [omega;0];%Check this for the Attitude Control of the Chaser
    psi = [p(1)*omega(2)*omega(3); p(2)*omega(1)*omega(3);p(3)*omega(1)*omega(2)];
    q_dot = 0.5*fn_CrossTensor(q_omega_rel,0)*q;
    tau = [tau_1;tau_2;tau_3];
    e_force = [e_force_x;e_force_y;e_force_z];
    omega_dot = psi + J_k*tau.*rand(3,1); %for static conditions
    %omega_dot = zeros(3,1);
    p_dot = zeros(3,1);
    rho_t_dot = zeros(3,1);
    ita_dot = zeros(4,1);
    r_ddot = fn_Compute_r_ddot(n,r,e_force,r_dot);
    %r_dot = zeros(3,1);
    %r_ddot = zeros(3,1);
    dy = [q_dot;omega_dot;p_dot;r_dot;r_ddot;rho_t_dot;ita_dot];
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
%Create state propagation matrix of the whole system

%Function: fn_Create_Phi()
%Inputs: X_a - 23-state vector,
%        n - angular velocity of spacecrafts,
%        t_delta - sampling time
%Outputs: [Phi] (State-propagation Matrix)
%Functionality: Generates the State Propagation matrix
%Author: Hrishik Mishra
function[Phi] = fn_Create_Phi(X_a,n,t_delta)
    Phi_t = fn_Create_Phi_t(n,t_delta);
    Phi_r = fn_Create_Phi_r(X_a);
    Phi_theta = eye(6,6);
    Phi = [Phi_r, zeros(9,6),zeros(9,6);zeros(6,9),Phi_t,zeros(6,6);zeros(6,9),zeros(6,6),Phi_theta];

end
%Create state propagation matrix of the whole system

%Function: fn_Create_Phi_r()
%Inputs: X_a - 23-state vector,
%Outputs: [Phi_r] (State-propagation Matrix for Rotation component)
%Functionality: Generates the State Propagation matrix for Rotation
%component
%Author: Hrishik Mishra
function [Phi_r] = fn_Create_Phi_r(X_a)
    global t_delta;
    p = X_a(8:10);
    omega = X_a(5:7);
    M = fn_Create_M(p,omega);
    N = fn_Create_N(omega);
    
    A = [-fn_VectorToSkewSymmetricTensor(omega),0.5*eye(3,3),zeros(3,3);zeros(3,3),M,N;zeros(3,3),zeros(3,3),zeros(3,3)];
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
%Function: fn_Create_N()
%Inputs: omega - angular velocity between Chaser Grasping frame and Chaser
%        CoM
%Outputs: N (Linearized model matrix for inertial parameters)
%Functionality: Generates M
%Author: Hrishik Mishra
function N = fn_Create_N(omega)
%
    N = [omega(2)*omega(3), 0, 0; 0, omega(1)*omega(3), 0; 0, 0, omega(1)*omega(2)];
end

%% Not being used for now

%Create the State Transition Matrix 11 for rotation kinematics
%Inputs: omega: angular velocity vector in current time step
%        t_delta: sampling time.
function Phi_r11 = fn_Create_phi_r11(omega, t_delta)
%
    Omega_Tensor = fn_VectorToSkewSymmetricTensor(omega);
    omega_norm = norm(omega);
    if (omega_norm ~= 0)
        Phi_r11 = eye(3,3) - (sin(omega_norm*t_delta)/omega_norm)*Omega_Tensor + ((1-cos(omega_norm*t_delta))/omega_norm^2)*Omega_Tensor^2;
    else
        Phi_r11 = eye(3,3);
    end
end

function Phi_r12  = fn_Create_phi_r12(omega,M,t_delta)
%
    Phi_r12 = zeros(3,3);
    omega_flag = 1;
    if ( norm(omega) == 0 )
        omega_flag = 0;
    end
    if ( omega_flag == 1 )       
        v_lambda = real(eig(M));
        %Set up matrix to solve for Gamma
        A = [1 v_lambda(1) v_lambda(1)^2;1 v_lambda(2) v_lambda(2)^2;1 v_lambda(3) v_lambda(3)^2];
        Gamma = real(inv(A));
        Omega_Tensor = fn_VectorToSkewSymmetricTensor(omega);
        omega_norm = norm(omega);
        for l=1:3
            for j=1:3
                for k= 1:3
                    Phi_r12 = Phi_r12 + Gamma(l,j)*fn_phi12_jk(k,v_lambda(j),t_delta, omega_norm)*(Omega_Tensor^(k-1))*(M^(l-1));
                end
            end
        end
        Phi_r12 = 0.5*Phi_r12;
    end
end   

function phi_jk = fn_phi12_jk(k,lambda,t_delta,omega_norm)
    switch(k)
        case 1
            phi_jk = fn_phi12_j1(lambda,t_delta);
        case 2
            phi_jk = fn_phi12_j2(lambda, omega_norm,t_delta);
        case 3
            phi_jk = fn_phi12_j3(lambda, omega_norm, t_delta);
        otherwise
            phi_jk = 0;
    end
end

%Following functions are for jth lambda
function phi_j1 = fn_phi12_j1(lambda,t_delta)    
    phi_j1 = (lambda^-1)*(exp(lambda*t_delta) - 1);
end
function phi_j2 = fn_phi12_j2(lambda, omega_norm, t_delta)
    phi_j2 = ((lambda^2*omega_norm + omega_norm^3)^-1)*(omega_norm*cos(omega_norm*t_delta) + lambda*sin(omega_norm*t_delta) - exp(lambda*t_delta)*omega_norm);
end
function phi_j3 = fn_phi12_j3(lambda, omega_norm, t_delta)
    phi_j3 = (omega_norm^-2)*(lambda^-1) + (lambda^3*omega_norm^2 + lambda*omega_norm^4)*((lambda^2)*cos(omega_norm*t_delta) - omega_norm*lambda*sin(omega_norm*t_delta) + (omega_norm^2)*exp(omega_norm*t_delta));
end

function Phi_r13 = fn_Create_phi_r13(omega, t_delta,M,N)
%
    Phi_r13 = zeros(3,3);
    v_lambda = real(eig(M));
    A = [1 v_lambda(1) v_lambda(1)^2;1 v_lambda(2) v_lambda(2)^2;1 v_lambda(3) v_lambda(3)^2];
    omega_flag = 1;
    if ( norm(omega) == 0 )
        omega_flag = 0;
    end
    if ( omega_flag == 1)
        Gamma = real(inv(A));
        for i = 1:3
            for j = 1:3
                Phi_r13 = Phi_r13 + Gamma(i,j)*(v_lambda(j)^-1)*(exp(v_lambda(j)*t_delta)-1)*(M^(i-1))*N;
            end
        end
    end
end

function Phi_r22 = fn_Create_phi_r22(omega,t_delta,M)
%
    Phi_r22 = eye(3,3);
    v_lambda = real(eig(M));
    A = [1 v_lambda(1) v_lambda(1)^2;1 v_lambda(2) v_lambda(2)^2;1 v_lambda(3) v_lambda(3)^2];
    omega_flag = 1;
    if ( norm(omega) == 0 )
        omega_flag = 0;
    end
    if ( omega_flag == 1 )
        Gamma = real(inv(A));
        for i = 1:3
            for j = 1:3
                Phi_r22 = Phi_r22 + Gamma(i,j)*exp(v_lambda(j)*t_delta)*M^(i-1);
            end
        end
    end
end
    
function Phi_r23 = fn_Create_phi_r23(omega, t_delta, M, N)
%
    Phi_r23 = zeros(3,3);
    omega_flag = 1;
    if (norm(omega) == 0)
        omega_flag = 0;
    end
    if ( omega_flag == 1 )
        Omega_Tensor = fn_VectorToSkewSymmetricTensor(omega);
        omega_norm = norm(omega);
        v_lambda = real(eig(M));
        A = [1 v_lambda(1) v_lambda(1)^2;1 v_lambda(2) v_lambda(2)^2;1 v_lambda(3) v_lambda(3)^2];
        Gamma = real(inv(A));
        for i= 1:3
            for j = 1:3
                for k = 1:3
                    Phi_r23 = Phi_r23 + Gamma(i,j)*fn_phi23_jk(k,v_lambda(j),omega_norm,t_delta)*(Omega_Tensor^(k-1))*(M^(i-1))*N;
                end
            end
        end
    end
end
function phi_jk = fn_phi23_jk(k,lambda,omega_norm,t_delta)
    switch(k)
        case 1
            phi_jk = fn_phi23_j1(lambda, t_delta);
        case 2
            phi_jk = fn_phi23_j2(lambda, omega_norm,t_delta);
        case 3
            phi_jk = fn_phi23_j3(lambda, omega_norm, t_delta);
        otherwise
            phi_jk = 0;
    end
end

function phi_j1 = fn_phi23_j1(lambda, t_delta)
    phi_j1 = (lambda^-2)*(1-exp(lambda*t_delta)) - (lambda^-1)*t_delta;
end

function phi_j2 = fn_phi23_j2(lambda, omega_norm, t_delta)
    phi_j2 = ((lambda^3)*(omega_norm^2) + (lambda*(omega_norm^4))^-1)*(omega_norm*lambda*sin(omega_norm*t_delta) + (lambda^2)*cos(omega_norm*t_delta) + (omega_norm^2)*exp(omega_norm*t_delta) - lambda^2 - omega_norm^2);
end
    
function phi_j3 = fn_phi23_j3(lambda, omega_norm, t_delta)
    phi_j3 = (omega_norm^-2)*(lambda^-2) + (((lambda^4)*(omega_norm^3) + (lambda^2)*(omega_norm^5))^-1)*((lambda^3)*sin(omega_norm*t_delta) - omega_norm*(lambda^2)*cos(omega_norm*t_delta) - (omega_norm^3)*exp(omega_norm*t_delta) - ((omega_norm^3)*lambda + (lambda^3)*omega_norm)*t_delta);
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

%Function: fn_Create_Q_r()
%Inputs: Phi - State Propagation Matrix for Rotational component
%        X - 23-state vector
%        t_delta - sampling time
%        tau - Torque disturbance
%        sig_p - Uncertainty in Inertial ratios (p)
%Outputs: Q_r (Rotation component Process noise)
%Functionality: Generates Q_r
%Author: Hrishik Mishra
function Q_r = fn_Create_Q_r(Phi,X, t_delta, tau,sig_p)
    p = X(7:9,1);
    J = fn_Create_J(p);
    Q_r11 = fn_Create_Q_r11(Phi(1:3,4:6), J, t_delta, tau);
    Q_r12 = fn_Create_Q_r12(Phi(1:3,4:6),Phi(4:6,4:6),J,t_delta,tau);
    Q_r22 = fn_Create_Q_r22(Phi(4:6,4:6),J,t_delta,tau);
    Q_p =  sig_p*eye(3,3);
    Q_r = [Q_r11, Q_r12, zeros(3,3);Q_r12' Q_r22 zeros(3,3);zeros(3,3) zeros(3,3) Q_p];
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

%Function: fn_Create_J()
%Inputs: p: Inertial ratios
%Outputs: J (Inertial Ratio Matrix for omega_dot)
%Functionality: Generates J
%Author: Hrishik Mishra
function J = fn_Create_J(p)
    J = [1 0 0;0 (1-p(2))/(1+p(1)) 0;0 0 (1+p(3))/(1-p(1))];
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

%Function: fn_Create_H()
%Inputs: q_nominal: Quaternion between Chaser CoM and Target CoM
%        X: 21-state vector
%Outputs: H (Linearized measurement function matrix)
%Functionality: Generates H
%Author: Hrishik Mishra
function H = fn_Create_H(q_nominal,X)
    rho_t_k = X(16:18,1);
    Matrix_1 = fn_Create_del_h1(q_nominal,rho_t_k);
    del_ita_v_k = X(19:21,1);
    del_q_v_k = X(1:3,1);
    Matrix_2 = fn_Create_del_h2(del_ita_v_k,del_q_v_k);
    H = [Matrix_1;Matrix_2];
    function Matrix_1 = fn_Create_del_h1(q,rho_t_k)
        q_0 = q(4);
        q_v = q(1:3);
        Q_v = fn_VectorToSkewSymmetricTensor(q_v);
        R = (2*q_0^2 - 1)*eye(3,3) + 2*q_0*Q_v + 2*q_v*q_v';

        first = -2*R*fn_VectorToSkewSymmetricTensor(rho_t_k);
        Matrix_1 = [first, zeros(3,6),eye(3,3),zeros(3,3),R,zeros(3,3)];
    end
    function Matrix_2 = fn_Create_del_h2(del_ita_v_k,del_q_v_k)    
        del_ita_0_k = sqrt(1 - norm(del_ita_v_k));
        del_q_0_k = sqrt(1 - norm(del_q_v_k));
        del_ita_k = [del_ita_v_k;del_ita_0_k];
        del_q_k = [del_q_v_k;del_q_0_k];
        del_Ita_k = fn_VectorToSkewSymmetricTensor(del_ita_v_k);
        first = -del_Ita_k + del_ita_k(4)*eye(3,3);
        del_Q_k = fn_VectorToSkewSymmetricTensor(del_q_v_k) ;
        last = del_Q_k + del_q_k(4)*eye(3,3);
        Matrix_2 = [first, zeros(3,6),zeros(3,3),zeros(3,3),zeros(3,3),last];
    end
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
    zk(4:6,1) = v(1:3);
end

%Function: fn_Create_h()
%Inputs: q: quaternion between Chaser CoM and Target CoM
%        X: 21-state vector 
%Outputs: h (non-linear measurement function)
%Functionality: Generates h
%Author: Hrishik Mishra
function h = fn_Create_h(q,X)
%
    r = X(10:12,1);
    rho_t = X(16:18,1);
    del_q_v = X(1:3,1);
    del_ita_v = X(19:21,1);
    del_ita_0 = sqrt(1 - norm(del_ita_v));
    del_q_0 = sqrt(1 - norm(del_q_v));
    del_ita = [del_ita_v;del_ita_0];
    del_q = [del_q_v;del_q_0];
    q_0 = q(4);
    q_v = q(1:3);
    Q_v = fn_VectorToSkewSymmetricTensor(q_v);
    R_q = (2*q_0^2 - 1)*eye(3,3) + 2*q_0*Q_v + 2*(q_v)*(q_v)';
    h1 = r + R_q*rho_t;
    temp = fn_CrossTensor(del_ita,0)*del_q;
    h2 = temp(1:3);
    h = [h1;h2];    
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