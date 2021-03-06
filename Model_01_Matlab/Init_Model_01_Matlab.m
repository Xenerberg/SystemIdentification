%Init the Model_01 (Matlab)
i_xx = 4; %Kg.m^2
i_yy = 8;
i_zz = 5;
InM = diag([i_xx,i_yy,i_zz]);
I_body = diag([i_xx;i_yy;i_zz]);

p_x = (i_yy - i_zz)/i_xx;
p_y = (i_zz - i_xx)/i_yy;
p_z = (i_xx - i_yy)/i_zz;


sig_tau = 2.5e-5; %rad^2/s^4
sig_p = 1e-5;
sig_theta = 1e-5;
sig_f = 2e-5;

%Initial conditions
%Quaternions
q_0_0 = 1;
q_1_0 = 0;
q_2_0 = 0;
q_3_0 = 0;
%omega
om_x_0 = 0.0;
om_y_0 = 0.0;
om_z_0 = 0.0;
%Inertial parameters
%Incorrect values
p_x_0 = 0;%2;
p_y_0 = 0;%-0.5;
p_z_0 = 0;%-2;
%Set correct values
%p_x_0 = p_x + 0.2*p_x;
%p_y_0 = p_y - 0.2*p_y;
%p_z_0 = p_z + 0.1*p_z;

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
%Set incorrect value
rho = [0;0;0];
%rho = [0.14;0.2;0.55];


%%Target point orientation
%ita =   [0.1029;0.1029; -0.2059; 0.9677];
ita = [0;0;0;1];

%Measurement data

Cov_r = 1e-1*eye(3,3);%3e-3*eye(3,3);
Cov_nu = 5e-1*eye(4,4);

%Orbital parameters
n = 0.0012;

%Simulation parameters
global t_delta; t_delta = .1;

X_a_0 = [q_1_0;q_2_0;q_3_0;q_0_0;om_x_0;om_y_0;om_z_0;p_x_0;p_y_0;p_z_0;r_x_0;r_y_0;r_z_0;r_dot_x_0;r_dot_y_0;r_dot_z_0;rho;ita];



rho_c = [0;0;0];

