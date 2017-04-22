close all;
figure;
title('q of body');
for iCount = 1:4
   subplot(2,2,iCount);
   plot(t,eval(strcat('q_',num2str(iCount-1))));
   ylabel((strcat('q_',num2str(iCount-1))));
end


figure;
title('omega in Inertial');
subplot(3,1,1);
plot(t,om_x);
ylabel('om_i_x');
subplot(3,1,2);
plot(t,om_y);
ylabel('om_i_y');
subplot(3,1,3);
plot(t,om_z);
ylabel('om_i_z');

figure;
title('rho_t in body');
subplot(3,1,1);
plot(t,rho_t_x);
ylabel('rho_t_x');
subplot(3,1,2);
plot(t,rho_t_y);
ylabel('rho_t_y');
subplot(3,1,3);
plot(t,rho_t_z);
ylabel('rho_t_z');

figure;
title('r in inertial');
subplot(3,1,1);
plot(t,r_x);
ylabel('r_x');
subplot(3,1,2);
plot(t,r_y);
ylabel('r_y');
subplot(3,1,3);
plot(t,r_z);
ylabel('r_z');

figure;
title('r_dot in inertial');
subplot(3,1,1);
plot(t,r_dot_x);
ylabel('r_dot_x');
subplot(3,1,2);
plot(t,r_dot_y);
ylabel('r_dot_y');
subplot(3,1,3);
plot(t,r_dot_z);
ylabel('r_dot_z');

figure;
title('omega in body');
subplot(3,1,1);
plot(t,om_b_x);
ylabel('om_b_x');
subplot(3,1,2);
plot(t,om_b_y);
ylabel('om_b_y');
subplot(3,1,3);
plot(t,om_b_z);
ylabel('om_b_z');

