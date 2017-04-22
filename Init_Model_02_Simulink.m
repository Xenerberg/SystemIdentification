%Init the Model_01 (Matlab)
i_xx = 4; %Kg.m^2
i_yy = 8;
i_zz = 5;

p_x = (i_yy - i_zz)/i_xx;
p_y = (i_zz - i_xx)/i_yy;
p_z = (i_xx - i_yy)/i_zz;
p = [p_x;p_y;p_z];

sig_tau = 2.5e-5; %rad^2/s^4
sig_p = 1e-4;
sig_theta = 1e-6;
sig_f = 2e-6;

%Initial conditions
%Quaternions
q_0_0 = 1;
q_1_0 = 0;
q_2_0 = 0;
q_3_0 = 0;
%omega
om_x_0 = 0;
om_y_0 = 0.0;
om_z_0 = 0;
%Inertial parameters
p_x_0 = 2;
p_y_0 = -0.5;
p_z_0 = -2;
%position
r_x_0 = 0;
r_y_0 = 0;
r_z_0 = 0;
%velocity
r_dot_x_0 = 0;
r_dot_y_0 = 0;
r_dot_z_0 = 0;


%Force perterubation
e_force_x = 0;
e_force_y = 0;
e_force_z = 0;

%Torque perturbation
tau_1 = 0;
tau_2 = 0;
tau_3 = 0;

%Target point position
rho = [0.2;0.1;0.05];
%Target point orientation
ita = [0.12;0.05;-0.15;0.98];
%ita = [0;0;0;1];

%Measurement data

Cov_r = 3e-3*eye(3,3);
Cov_nu = 5e-3*eye(4,4);



%Simulation parameters
global t_delta; t_delta = 0.1;
global t_Kalman; t_Kalman = 0.01;
global t_Sim; t_Sim = 0.001;
global totalSimulationTime;totalSimulationTime = 200;
%Orbital parameters
global parameter_gravitation;parameter_gravitation = 398.6005e12;
global n;n = 0.0012; %Orbital velocity of Chaser/Target (nearly same)

X_a_0 = [q_1_0;q_2_0;q_3_0;q_0_0;om_x_0;om_y_0;om_z_0;r_x_0;r_y_0;r_z_0;r_dot_x_0;r_dot_y_0;r_dot_z_0];

q = struct('q_1',q_1_0,'q_2',q_2_0,'q_3',q_3_0,'q_0',q_0_0);
om = struct('omega_x',om_x_0,'omega_y',om_y_0,'omega_z',om_z_0);
r = struct('r_x',r_x_0,'r_y',r_y_0,'r_z',r_z_0);
r_dot = struct('r_x_dot',r_dot_x_0,'r_y_dot',r_dot_y_0,'r_z_dot',r_dot_z_0);
X_a_0_struct = struct('q',q,'omega',om,'r',r,'r_dot',r_dot);
X_a_0_bus = Simulink.Bus.createObject(X_a_0_struct);

rho_c = [0;0;0.9];

