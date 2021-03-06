
%Program: DLR-project 
%Topic: Motion and Parameter estimation of space objects using Laser-Vision
%data (Extended Kalman filter/Quaternion algebra)
%Author: Hrishik Mishra

%Function: fn_Main()
%Inputs: none
%Outputs:none
%Functionality: Initialize the system/Execute the Kalman filter/Display results
%Author: Hrishik Mishra

function [ ] = fn_Main( )
    %Implementation of the Motion parameter estimation
    clc; %close all;;close all;
    %Global variables
    global n;n = 0.0012; %Orbital velocity of Chaser/Target (nearly same)
    global parameter_gravitation;parameter_gravitation = 398.6005e12;
    global totalSimulationTime;totalSimulationTime = 50;        
   
    
    %1) Initialize the system with model parameters
    Init_Model_02_Matlab;   
    
    %Integrate Dynamic equations to see result;    (Simulated data)
    %%
    dy_time = (0:t_delta:totalSimulationTime)';
    s1 = 0.05;
    s2 = 0.1;tdelay = 0.09;
    Time = dy_time;
    options = odeset('RelTol',1e-7,'AbsTol',1e-9*ones(13,1));
    X_a_0_sim = load('./Model_01_Matlab/X_a_0_Model_02');
    [T,X_a] = ode45(@fn_StateSpace,dy_time,X_a_0_sim.X_a_0_Model_02,options);
    Signal_Quaternion = [];
    Mu = (fn_CrossTensor(ita,0)*X_a(:,1:4)')';
    Mu_noise = zeros(size(Mu));
    Mu_noise = [Mu_noise(:,4), Mu_noise(:,1:3)];
    Mu_noise = quatnormalize(Mu + normrnd(0,sqrt(1e-3),size(Mu)));
    Mu_noise = [Mu_noise(:,1:3), Mu_noise(:,4)];
    Signal_Quaternion(:,2:4) = Mu_noise(:,1:3);
    Signal_Quaternion(:,1) = Mu_noise(:,4); %+ %0.01*randn(length(X_a),1);%0.01*randn(length(X_a),1);
    
    sigma_z = 3e-5;
    r_c = zeros(3,length(dy_time));    
    Mu_noise = Mu_noise';
    r_c = zeros(3,length(dy_time));
    for iCount = 1:length(X_a)
        r_c(1:3,iCount) = X_a(iCount,8:10)' + rho_c + fn_CreateRotationMatrix(X_a(iCount,1:4)')*rho;
    end
    r_c_noise(1:3,:) = r_c(1:3,:) + normrnd(0,sqrt(1e-3),size(r_c));
    Signal_vector = timeseries([r_c_noise(1:3,:);Signal_Quaternion'],T,'Name','Signal');
    meas_1 = resample(Signal_vector,1:s1:totalSimulationTime);
    Mu_noise = [Mu_noise(:,4), Mu_noise(:,1:3)];
    Mu_noise = quatnormalize(Mu + normrnd(0,sqrt(1e-5),size(Mu)));
    Mu_noise = [Mu_noise(:,1:3), Mu_noise(:,4)];
    %Mu_noise = circshift(Mu_noise,[tdelay,0]);
    Signal_Quaternion_(:,2:4) = Mu_noise(:,1:3);
    Signal_Quaternion_(:,1) = Mu_noise(:,4); %+ %0.01*randn(length(X_a),1);%0.01*randn(length(X_a),1)
    for iCount = 1:length(X_a)
        r_c(1:3,iCount) = X_a(iCount,8:10)' + rho_c + fn_CreateRotationMatrix(X_a(iCount,1:4)')*rho;
    end
    r_c_noise(1:3,:) = r_c(1:3,:) + normrnd(0,sqrt(1e-5),size(r_c));
    %r_c_noise = circshift(r_c_noise,[0,tdelay]);
    Signal_vector_ = timeseries([r_c_noise(1:3,:);Signal_Quaternion_'],T,'Name','Signal');
    s2T = 1:s2:totalSimulationTime;
    meas_2 = resample(Signal_vector_,s2T);    
    meas_2.Time = meas_2.Time + tdelay;
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
    X_a_init = X_a_0; %23-state vector to store the final output of Integration of dynamic equations (Predicted Result)
    X_a_pre = zeros(13, length(Time));%Time-series of X_a_init (23-state vector)
    X_a_itr = X_a_0;%23-state vector to store the output after Update-Stage (Estimated result)    
    X_a_Estimated = zeros(13, length(Time));%Time-series of X_a_itr
    X_a_Estimated(:,1) = X_a_itr;%Initialize the vector
    X_pre = zeros(12,1);%21-state vector to store the output after Prediction-Stage
    X_pre_estimate = zeros(12,1);%Time-series of 21-state vector after Update-stage    
    X_a_pre(:,1) = X_a_init;
    v_time = dy_time; %Time vector   
    X_pre(4:12,1) = X_a_0(5:13,1);
    tr_P = zeros(1,length(Time));
    %Describe the state-error covariance initial matrix
    P_post = eye(12,12);    
    C_pre = eye(12,12);
    C_post = eye(12,12);
    tr_C(1) = trace(C_post);
    delta_conv = zeros(1,length(X_a_Estimated));
    delta_conv = 0;
%     P_post(1:3,1:3) = 0.5*eye(3,3);  
%     P_post(7:9,7:9) = 4*eye(3,3);
%     P_post(10:12,10:12) = eye(3,3);    
    rankObs = zeros(length(v_time));
    residuals = zeros(length(v_time)-1,6);
    sq_error = zeros(length(v_time),1);
    error_post = X_a(iCount-1,:)' -  X_a_Estimated(:,iCount-1);
    sq_error(1) = error_post'*error_post;
    tr_M(1) = trace(error_post*error_post');
    post_residual = residuals;
    signal = zeros(7,1);
    test = zeros(6,length(v_time));
    Mh_PostDist = zeros(length(v_time),1);
    Mh_PreDist = Mh_PostDist;
    %%
    q_nominal_1 = X_a_itr(1:4);
    q_nominal_2 = X_a_itr(1:4);
    meas_1count = 1;
    nexttime_1 = meas_1.Time(meas_1count);
    meas_2count = 1;
    nexttime_2 = meas_2.Time(meas_2count);
    nextstime_2 = s2T(meas_2count);
    oosmFlag = false;initFlag = true;
    %%Kalman filter implementation in an iterative loop
    tic;
    for iCount = 2:length(v_time)
        if abs(dy_time(iCount) - nextstime_2) < 1e-4 
            %oosmFlag = true;                 
        end
        %q_nominal = nominal-quaternion (Chaser CoM-Target CoM) from previous states
                
        
        %% Run the dynamic model to get predicted states        
        [X] = fn_StateSpace(0,X_a_Estimated(:,iCount-1));
        X_a_init = X_a_Estimated(:,iCount-1) + t_delta*X;        
        X_a_pre(:,iCount) = X_a_init;        
        X_a_Estimated(:,iCount) = X_a_pre(:,iCount);
%         [~,X] = ode45(@fn_StateSpace,[v_time(iCount-1),v_time(iCount)],X_a_Estimated(:,iCount-1),options);                
%         X_a_pre = X(end,:)';
%         X_a_Estimated(:,iCount) = X_a_pre;
        
        
        %% State propagation               
        A = fn_Create_A(X_a_itr,n,p);
        Phi = fn_Create_Phi(X_a_itr, n, t_delta,p);
        Phi_C = fn_Create_Phi(X_a(iCount,:)', n, t_delta,p);
        
        %% Create Process-noise covariance matrix
        Q_r = fn_Create_Q_r(Phi,X_a_itr, t_delta, sig_tau,sig_p,p);
        Q_t  = fn_Create_Q_t_k(t_delta, n,sig_f);
        Q = [Q_r, zeros(6,6);zeros(6,6),Q_t];
        %Q = 0.1*eye(12,12);
        
        %% Create State-error Covariance for Prediction-Stage
        P_pre = Phi*P_post*Phi' + Q;
        C_pre = Phi_C*C_post*Phi_C' + Q;
        
        P_post = P_pre;
        C_post = C_pre;
        K = zeros(12,6);
        
        zk = zeros(6,1);
        h = zeros(6,1);
        alpha = 0;
        if abs(dy_time(iCount) - nexttime_1) < 1e-4
            %% Change to intermediate state
            del_q = fn_CrossTensor(X_a_init,0)*[-q_nominal_1(1:3);q_nominal_1(4)];  
            %fprintf('[%f %f %f %f]\n',del_q(1),del_q(2),del_q(3),del_q(4));
            %Set up X_pre for the Update-stage
            X_pre(1:3) = del_q(1:3); 
            X_pre(4:12) = X_a_init(5:13);

            %% Fetch measurement
            signal(1:3,1) = meas_1.Data(1:3,:,meas_1count); 
            signal(4:6,1) = meas_1.Data(5:7,:,meas_1count); 
            signal(7,1) = meas_1.Data(4,:,meas_1count);
%             signal(1:3,1) = Signal_vector.Data(1:3,:,iCount); 
%             signal(4:6,1) = Signal_vector.Data(5:7,:,iCount); 
             
            if meas_1.Time(meas_1count) < totalSimulationTime
                meas_1count = meas_1count+1;
                nexttime_1 = meas_1.Time(meas_1count);
            end
            %signal(4:7,1) = [q_measured(2:4); q_measured(1)];
            zk =fn_Create_obs(signal,rho_c,q_nominal_1,ita);
            fprintf('[%f %f %f]\n',zk(4),zk(5),zk(6));

            %% Update phase
            Hk = fn_Create_H(q_nominal_1,rho);
            Sk = fn_Create_S(q_nominal_1,ita,Cov_r_1,Cov_nu_2); %+ 5*Hk*P_pre*Hk';
            %Bump up EKF        
            %Sk_new = Sk + Hk*P_pre*Hk';
            %K = fn_ComputeKalmanGain(P_pre,Hk,Sk_new);
            K = fn_ComputeKalmanGain(P_pre,Hk,Sk);
            S = Hk*P_pre*Hk' + Sk;
            h = fn_Create_h(fn_CrossTensor(del_q,0)*q_nominal_1,X_pre, rho); 
            %h = fn_Create_h(X_a_init(1:4),X_pre, rho);        
            residuals(meas_1count-1,:) = zk - h;
            
            Mh_PreDist(meas_1count-1) = sqrt(residuals(meas_1count-1,:)*inv(Hk*P_pre*Hk' + Sk)*residuals(meas_1count-1,:)');
            X_pre_estimate = X_pre(:,end) + K*(zk - h);       


            %% Error/Warning checks in computed Quaternions
            del_q_v = X_pre_estimate(1:3);
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



            %% Convert from intermediate to Estimated-states after Update-stage 
            q_est = fn_CrossTensor(del_q,0)*q_nominal_1;              
            X_a_itr = [q_est;X_pre_estimate(4:12)];
            X_a_Estimated(:,iCount) = X_a_itr;
            P_post = (eye(12,12)-K*Hk)*P_pre;
            
            post_h = fn_Create_h(X_a_Estimated(1:4,iCount),X_pre_estimate(:),rho);
            post_residual(meas_1count-1,:) = zk - post_h;
            Hk_post = fn_Create_H(q_nominal_1,rho);
            Mh_PostDist(meas_1count-1) = sqrt(post_residual(meas_1count-1,:)*inv(Hk_post*P_post*Hk_post' + Sk)*post_residual(meas_1count-1,:)');
            %[eig_vectorObs,eig_Observer] = eig(A - K*Hk);
            

            Hk_true = fn_Create_H(X_a(iCount,1:4),rho);
            Sk_true = sigma_z*eye(6,6);
            C_post = (eye(12,12) - fn_ComputeKalmanGain(C_pre,Hk_true,Sk_true)*Hk_true)*C_pre;
            %C_pre - (C_pre*Hk_true'/(Hk_true*C_pre*Hk_true' + Sk))*Hk_true*C_pre' + 1e-1*eye(12,12);
            tr_C(meas_1count) = trace(C_post);
            delta_conv(meas_1count) =  max(eig((Sk)))*max(eig(inv(Sk)*Hk*P_pre*Hk'*inv(Hk*P_pre*Hk'+Sk)));
            eigMin(meas_1count) = min(real(eig(A-K*Hk)));
            
            q_nominal_1 = expm(0.5*fn_CrossTensor([X_a_Estimated(5:7,iCount);0],0)*s1)*X_a_Estimated(1:4,iCount);
            q_nominal = q_nominal_1;
            alpha = 1;
        end
        
        %q_nominal_2 = q_nominal_1;
        %P_pre = P_post;        
        C_pre = C_post;
        
        
        %OOSM
        if oosmFlag == true
            if initFlag == true
                initFlag = false;
                del_q_oosm = fn_CrossTensor(X_a_init(1:4),0)*[-q_nominal(1:3);q_nominal(4)];  
                %fprintf('[%f %f %f %f]\n',del_q(1),del_q(2),del_q(3),del_q(4));
                %Set up X_pre for the Update-stage
                X_pre_oosm(1:3) = del_q_oosm(1:3); 
                X_pre_oosm(4:12) = X_a_init(5:13);   
                U = (eye(12,12) - K*Hk)*Phi*P_pre;
                P_d = P_pre - alpha*P_pre*Phi'*Hk'*inv(S)*Hk*Phi*P_pre;
            else                
                x_d = X_pre_oosm' + U'*Phi'*Hk'*inv(S)*(zk-h);
                P_d = P_d - alpha*U'*Phi'*Hk'*inv(S)*Hk*Phi*U;
                U = (eye(12,12) - K*Hk)*Phi*U;                
            end
        end

       
        if abs(dy_time(iCount) - nexttime_2) < 1e-4
            %% Change to intermediate state
            del_q = fn_CrossTensor(X_a_init(1:4),0)*[-q_nominal_2(1:3);q_nominal_2(4)];  
            %fprintf('[%f %f %f %f]\n',del_q(1),del_q(2),del_q(3),del_q(4));
            %Set up X_pre for the Update-stage
            X_pre(1:3) = del_q(1:3); 
            X_pre(4:12) =X_a_init(5:13);

            %% Fetch measurement
            signal(1:3,1) = meas_2.Data(1:3,:,meas_2count); 
            signal(4:6,1) = meas_2.Data(5:7,:,meas_2count); 
            signal(7,1) = meas_2.Data(4,:,meas_2count);
%             signal(1:3,1) = Signal_vector.Data(1:3,:,iCount); 
%             signal(4:6,1) = Signal_vector.Data(5:7,:,iCount); 
             
            if meas_2.Time(meas_2count) < totalSimulationTime
                meas_2count = meas_2count+1;
                nexttime_2 = meas_2.Time(meas_2count);
                nextstime_2 = s2T(meas_2count);
            end
            %signal(4:7,1) = [q_measured(2:4); q_measured(1)];
            zk =fn_Create_obs(signal,rho_c,q_nominal_2,ita);
            fprintf('[%f %f %f]\n',zk(4),zk(5),zk(6));

            %% Update phase
            Hk = fn_Create_H(q_nominal_2,rho);
            Sk = fn_Create_S(q_nominal_2,ita,Cov_r_2,Cov_nu_2); %+ 5*Hk*P_pre*Hk';
            %Bump up EKF        
            %Sk_new = Sk + Hk*P_pre*Hk';
            %K = fn_ComputeKalmanGain(P_pre,Hk,Sk_new);
            K = fn_ComputeKalmanGain(P_pre,Hk,Sk);
            %K = fn_ComputeKalmanGain(U,Hk,Sk);
            
            h = fn_Create_h(fn_CrossTensor(del_q,0)*q_nominal_2,X_pre, rho); 
            %h = fn_Create_h(X_a_init(1:4),X_pre, rho);        
            residuals(meas_2count-1,:) = zk - h;
            
            Mh_PreDist(meas_2count-1) = sqrt(residuals(meas_2count-1,:)*inv(Hk*P_pre*Hk' + Sk)*residuals(meas_2count-1,:)');
            X_pre_estimate = X_pre(:,end) + K*(zk - h);       

            
            %% Error/Warning checks in computed Quaternions
            del_q_v = X_pre_estimate(1:3);
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



            %% Convert from intermediate to Estimated-states after Update-stage 
            q_est = fn_CrossTensor(del_q,0)*q_nominal_2;              
            X_a_itr = [q_est;X_pre_estimate(4:12)];
            X_a_Estimated(:,iCount) = X_a_itr;
            P_post = (eye(12,12)-K*Hk)*P_pre;
            %P_post = P_pre - K*inv(Hk*P_pre*Hk' + Sk)*K';
            
            post_h = fn_Create_h(X_a_Estimated(1:4,iCount),X_pre_estimate(:),rho);
            post_residual(meas_2count-1,:) = zk - post_h;
            Hk_post = fn_Create_H(q_nominal_2,rho);
            Mh_PostDist(meas_2count-1) = sqrt(post_residual(meas_2count-1,:)*inv(Hk_post*P_post*Hk_post' + Sk)*post_residual(meas_2count-1,:)');
            %[eig_vectorObs,eig_Observer] = eig(A - K*Hk);

            Hk_true = fn_Create_H(X_a(iCount,1:4),rho);
            Sk_true = sigma_z*eye(6,6);
            C_post = (eye(12,12) - fn_ComputeKalmanGain(C_pre,Hk_true,Sk_true)*Hk_true)*C_pre;
            %C_pre - (C_pre*Hk_true'/(Hk_true*C_pre*Hk_true' + Sk))*Hk_true*C_pre' + 1e-1*eye(12,12);
            
            delta_conv(meas_2count) =  max(eig((Sk)))*max(eig(inv(Sk)*Hk*P_pre*Hk'*inv(Hk*P_pre*Hk'+Sk)));
            eigMin(meas_2count) = min(real(eig(A-K*Hk)));
            
            q_nominal_2 = expm(0.5*fn_CrossTensor([X_a_Estimated(5:7,iCount);0],0)*s2)*X_a_Estimated(1:4,iCount);
            q_nominal = q_nominal_2;
        end
        tr_P(iCount) = trace(P_post);
        error_post = X_a(iCount,:)' -  X_a_Estimated(:,iCount);
        sq_error(iCount) = sqrt(error_post'*error_post);
        tr_M(iCount) = trace(error_post*error_post');
        tr_C(iCount) = trace(C_post);
        
    end
    toc;
    save('Data\trPi_TD_delS.mat','tr_M');
    %% Generate measurements to compare with actual measurements
    Mu_est = zeros(4,length(X_a_Estimated));
    r_c_est = zeros(3,length(X_a_Estimated));
    
    for iCount = 1:length(X_a_Estimated)
       Mu_est(:,iCount) = fn_CrossTensor(ita,0)*X_a_Estimated(1:4,iCount);
       r_c_est(1:3,iCount) = X_a_Estimated(8:10,iCount) + rho_c + fn_CreateRotationMatrix(X_a_Estimated(1:4,iCount))*rho;
       err_r_c(iCount) = norm(r_c_est(1:3,iCount) - r_c(1:3,iCount));
       temp = fn_CrossTensor(Mu(iCount,:)',0)*[-Mu_est(1:3,iCount);Mu_est(4,iCount)];
       err_Mu(iCount) = asind(norm(temp(1:3)));
    end
    
    %% Plot the details
    figure('units','normalized','outerposition',[0 0 1 1])
    subplot(4,2,1);
    stairs(dy_time,reshape(Signal_vector.Data(1,:,:),length(Signal_vector.Data(1,:,:)),1),'-');
    hold all;
    plot(dy_time,r_c_est(1,:),'LineWidth',2,'Color','red');
    legend('Measured','Estimated');
    ylabel('$r_{c_x}$','interpreter','latex','FontSize', 20);
    xlabel('$(a)$','interpreter','latex','FontSize', 15);
    subplot(4,2,2);
    stairs(dy_time,reshape(Signal_vector.Data(2,:,:),length(Signal_vector.Data(2,:,:)),1),'-');
    hold all;
    plot(dy_time,r_c_est(2,:),'LineWidth',2,'Color','red');
    legend('Measured','Estimated');
    ylabel('$r_{c_y}$','interpreter','latex','FontSize', 20);
    xlabel('$(b)$','interpreter','latex','FontSize', 15);
    subplot(4,2,3);
    stairs(dy_time,reshape(Signal_vector.Data(3,:,:),length(Signal_vector.Data(3,:,:)),1),'-');
    hold all;
    plot(dy_time,r_c_est(3,:),'LineWidth',2,'Color','red');
    legend('Measured','Estimated');
    ylabel('$r_{c_z}$','interpreter','latex','FontSize', 20);
    xlabel('$(c)$','interpreter','latex','FontSize', 15);
    subplot(4,2,5);
    stairs(dy_time,Signal_Quaternion(:,1));
    hold all;
    plot(dy_time,Mu_est(4,:),'LineWidth',2,'Color','red');
    legend('Measured','Estimated');
    ylabel('$\mu_0$','interpreter','latex','FontSize', 20);
    xlabel('$(d)$','interpreter','latex','FontSize', 15);
    subplot(4,2,6);
    stairs(dy_time,Signal_Quaternion(:,2));
    hold all;
    plot(dy_time,Mu_est(1,:),'LineWidth',2,'Color','red');
    legend('Measured','Estimated');
    ylabel('$\mu_1$','interpreter','latex','FontSize', 20);
    xlabel('$(e)$','interpreter','latex','FontSize', 15);
    subplot(4,2,7);
    stairs(dy_time,Signal_Quaternion(:,3));
    hold all;
    plot(dy_time,Mu_est(2,:),'LineWidth',2,'Color','red');
    legend('Measured','Estimated');
    ylabel('$\mu_2$','interpreter','latex','FontSize', 20);
    xlabel('$(f)$','interpreter','latex','FontSize', 15);
    subplot(4,2,8);
    stairs(dy_time,Signal_Quaternion(:,4))
    hold all;
    plot(dy_time,Mu_est(3,:),'LineWidth',2,'Color','red');
    legend('Measured','Estimated');
    ylabel('$\mu_3$','interpreter','latex','FontSize', 20);
    xlabel('$(g)$','interpreter','latex','FontSize', 15);
    print('Images\1','-depsc');
    
    
    figure('units','normalized','outerposition',[0 0 1 1])
    subplot(2,1,1);
    plot(dy_time,X_a_Estimated(5:7,:)');hold all;
    plot(dy_time,X_a(:,5:7),'LineStyle','-.');
    l = legend('$\hat{\omega_x}$','$\hat{\omega_y}$','$\hat{\omega_z}$','$\omega_x$','$\omega_y$','$\omega_z$');
    set(l,'Interpreter','latex','FontSize', 15);
    ylabel('$\omega$','interpreter','latex','FontSize', 20);
    xlabel('$(a)$','interpreter','latex','FontSize', 15);
    ylim([-4,4]);
    
%     subplot(3,1,2);
%     plot(dy_time,X_a_Estimated(8:10,:)');hold all;  
%     plot(dy_time,X_a(:,8:10),'LineStyle','-.')
%     legend('x_est','y_est','z_est','x_true','y_true','z_true');
%     ylabel('p');
    subplot(2,1,2);
    plot(dy_time,X_a_Estimated(11:13,:)');hold all;
    plot(dy_time,X_a(:,11:13),'LineStyle','-.');
    l = legend('$\hat{\dot{r}}_x$','$\hat{\dot{r}}_y$','$\hat{\dot{r}}_z$','$\dot{r}_x$','$\dot{r}_y$','$\dot{r}_z$');
    set(l,'Interpreter','latex','FontSize', 15);
    ylabel('$\dot{r}$','interpreter','latex','FontSize', 20);    
    xlabel('$(b)$','interpreter','latex','FontSize', 15);
    ylim([-0.05,0.05]);
    print('Images\2','-depsc');
    figure('units','normalized','outerposition',[0 0 1 1])
    subplot(3,2,1);
    plot(dy_time,X_a_Estimated(1,:)');hold all;
    plot(dy_time,X_a(:,1),'LineStyle','-.');
    xlabel('$(a)$','interpreter','latex','FontSize', 15);
    l = legend('$\hat{q_1}$','$q_1$');
    set(l,'Interpreter','latex','FontSize', 15);
    %ylabel('quaternion');
    subplot(3,2,2);    
    plot(dy_time,X_a_Estimated(2,:)');hold all;
    plot(dy_time,X_a(:,2),'LineStyle','-.');
    xlabel('$(b)$','interpreter','latex','FontSize', 15);
    l = legend('$\hat{q_2}$','$q_2$');
    set(l,'Interpreter','latex','FontSize', 15);
    
    subplot(3,2,3);
    plot(dy_time,X_a_Estimated(3,:)');hold all;
    plot(dy_time,X_a(:,3),'LineStyle','-.');
    xlabel('$(c)$','interpreter','latex','FontSize', 15);
    l = legend('$\hat{q_3}$','$q_3$');
    set(l,'Interpreter','latex','FontSize', 15);
    
    subplot(3,2,4);
    
    plot(dy_time,X_a_Estimated(4,:)');hold all;
    plot(dy_time,X_a(:,4),'LineStyle','-.');
    xlabel('$(d)$','interpreter','latex','FontSize', 15);
    l = legend('$\hat{q_0}$','$q_0$');
    set(l,'Interpreter','latex','FontSize', 15);
    
    subplot(3,2,5);
    plot(dy_time,X_a_Estimated(8:10,:));hold all;
    plot(dy_time,X_a(:,8:10),'LineStyle','-.');
    xlabel('$(e)$','interpreter','latex','FontSize', 15);
    l = legend('$\hat{r}_x$','$\hat{r}_y$', '$\hat{r}_z$', '$r_x$','$r_y$', '$r_z$');
    set(l,'Interpreter','latex','FontSize', 15);
    print('Images\3','-depsc');
    figure('units','normalized','outerposition',[0 0 1 1])
    subplot(4,1,1);
%     plot(dy_time,Mh_PreDist);
%     ylabel('$d[e_k(k+1|k)]$','interpreter','latex','FontSize', 15);
%     xlabel('$(a)$','interpreter','latex','FontSize', 15);
%     subplot(3,2,2);
%     plot(dy_time,Mh_PostDist);
%     ylabel('$d[e_k(k|k)]$','interpreter','latex','FontSize', 15);
%      xlabel('$(b)$','interpreter','latex','FontSize', 15);
%     subplot(3,2,3);
    plot(dy_time,(tr_P));hold all;
    plot(dy_time, tr_C);
    plot(dy_time, tr_M);
    ylabel('$trace$','interpreter','latex','FontSize', 15);
     xlabel('$(c)$','interpreter','latex','FontSize', 15);
    l = legend('$\Sigma(k|k)$','$CRB(k|k)$', '$\Pi(k|k)$');
    set(l,'Interpreter','latex','FontSize', 15);
    ylim([0,6e-2]);
    subplot(4,1,2)
    plot(dy_time, (sq_error));hold all;
    plot(dy_time, mean(sq_error)*ones(size(dy_time)),'LineStyle','-.');
    ylim([0,2e-1]);
    ylabel('$\tilde{x}(k|k)^T\tilde{x}(k|k)$','interpreter','latex','FontSize', 15);
    xlabel('$(d)$','interpreter','latex','FontSize', 15);
    
    subplot(4,1,3);    
    plot(dy_time,err_r_c);hold all;
    plot(dy_time, mean(err_r_c(200:end))*ones(size(dy_time)),'LineStyle','-.');
    ylabel('$e[m]$','interpreter','latex','FontSize', 15);
    xlabel('$(e)$','interpreter','latex','FontSize', 15);
    ylim([0,0.2]);
    subplot(4,1,4);
    plot(dy_time,err_Mu);hold all;
    plot(dy_time, mean(err_Mu(200:end))*ones(size(dy_time)),'LineStyle','-.');
    ylabel('$e[deg]$','interpreter','latex','FontSize', 15);
    xlabel('$(f)$','interpreter','latex','FontSize', 15);
    ylim([0,5]);
    print('Images\4','-depsc');
    %save('Model_01_Matlab\X_a_Estimated_1.mat','X_a_Estimated');
    %error = X_a' - X_a_Estimated;
    

%     figure;
%     %subplot(1,1,1);
%     plot(dy_time,eigMin);
%     ylabel('$\lambda_{min}$','interpreter','latex','FontSize', 15);
%     xlabel('$t(sec)$','interpreter','latex','FontSize', 15);
end

%Function: fn_StateSpace()
%Inputs: X_a (23-state vector)
%Outputs: [dy] (Time-derivative of the state vector)
%Functionality: Implements the State-space model of the dynamic system
%Author: Hrishik Mishra
function[dy] = fn_StateSpace(~,X_a)
    global n;
    Init_Model_02_Matlab;
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
%Create state propagation matrix of the whole system

%Function: fn_Create_Phi()
%Inputs: X_a - 23-state vector,
%        n - angular velocity of spacecrafts,
%        t_delta - sampling time
%Outputs: [Phi] (State-propagation Matrix)
%Functionality: Generates the State Propagation matrix
%Author: Hrishik Mishra
function[Phi] = fn_Create_Phi(X_a,n,t_delta,p)
    Phi_t = fn_Create_Phi_t(n,t_delta);
    Phi_r = fn_Create_Phi_r(X_a,p);
    Phi = [Phi_r, zeros(6,6);zeros(6,6),Phi_t];

end
%Create state propagation matrix of the whole system

%Function: fn_Create_Phi_r()
%Inputs: X_a - 23-state vector,
%Outputs: [Phi_r] (State-propagation Matrix for Rotation component)
%Functionality: Generates the State Propagation matrix for Rotation
%component
%Author: Hrishik Mishra
function [Phi_r] = fn_Create_Phi_r(X_a,p)
    global t_delta;    
    omega = X_a(5:7);
    M = fn_Create_M(p,omega);   
    
    A = [-fn_VectorToSkewSymmetricTensor(omega),0.5*eye(3,3);zeros(3,3),M];
    Phi_r = expm(A*t_delta);
%     phi_r11 = fn_Create_phi_r11(omega,t_delta);
%     phi_r12 = fn_Create_phi_r12(omega,M,t_delta);
%     phi_r22 = fn_Create_phi_r22(omega,t_delta,M);
%     phi_r13 = fn_Create_phi_r13(omega,t_delta,M,N);
%     phi_r23 = fn_Create_phi_r23(omega,t_delta,M,N);
%     Phi_r = [phi_r11, phi_r12, phi_r13;zeros(3,3),phi_r22,phi_r23;zeros(3,3),zeros(3,3),eye(3,3)];
%     

end

function [A] = fn_Create_A(X_a,n,p)
    omega = X_a(5:7);
    M = fn_Create_M(p,omega);   
    A_r = [-fn_VectorToSkewSymmetricTensor(omega),0.5*eye(3,3);zeros(3,3),M];
    K = [3*n^2 0 0;0 0 0;0 0 -n^2];
    A_t = [zeros(3,3),eye(3,3);K, -2*fn_VectorToSkewSymmetricTensor([0,0,n])];
    A = [A_r,zeros(6,6);zeros(6,6),A_t];

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
function h = fn_Create_h(q,X,rho)
%
    r = X(7:9,1);
    rho_t = rho;
    del_q_v = X(1:3,1);        
    del_q_0 = sqrt(1 - norm(del_q_v));    
    del_q = [del_q_v;del_q_0];
    q_0 = q(4);
    q_v = q(1:3);
    Q_v = fn_VectorToSkewSymmetricTensor(q_v);
    R_q = (2*q_0^2 - 1)*eye(3,3) + 2*q_0*Q_v + 2*(q_v)*(q_v)';
    h1 = r + R_q*rho_t;
    temp = del_q;
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