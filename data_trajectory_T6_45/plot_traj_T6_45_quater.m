%close all;
title('q of body');
for iCount = 1:4
   subplot(2,2,iCount);
   plot(time,eval(strcat('q_',num2str(iCount-1))));
   ylabel((strcat('q_',num2str(iCount-1))));
end

figure;
title('r in inertial');
subplot(3,1,1);
plot(time,r_x);
ylabel('r_x');
subplot(3,1,2);
plot(time,r_y);
ylabel('r_y');
subplot(3,1,3);
plot(time,r_z);
ylabel('r_z');