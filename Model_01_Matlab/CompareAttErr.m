%Script to compare attitude error representations
close all;
str_X_a_true = load('Model_01_Matlab/X_a_true.mat');
str_X_a_Est1 = load('Model_01_Matlab/X_a_Estimated_1.mat');
str_X_a_Est2 = load('Model_01_Matlab/X_a_Estimated_2.mat');

X_a_true = (str_X_a_true.X_a)';
X_a_Est1 = str_X_a_Est1.X_a_Estimated;
X_a_Est2 = str_X_a_Est2.X_a_Estimated;


figure;
title('quaternion details');
for iCount = 1:4
    
   subplot(2,2,iCount);
   plot(X_a_true(iCount,:),'linewidth',2,'linestyle','-');hold all;
   plot(X_a_Est1(iCount,:),'linewidth',2,'linestyle','--');
   plot(X_a_Est2(iCount,:),'linewidth',2,'linestyle',':');   
   xlim([0,500]);
   l = legend('true','small angle','gibbs vector');   
   switch(iCount)
       case 1
           ylabel('$q_1$','interpreter','latex','FontSize', 15);
       case 2
           ylabel('$q_2$','interpreter','latex','FontSize', 15);
       case 3
           ylabel('$q_3$','interpreter','latex','FontSize', 15);
       case 4
           ylabel('$q_0$','interpreter','latex','FontSize', 15);
   end
end



for iCount = 1:500
   temp1 = fn_CrossTensor(X_a_true(1:4,iCount),0)*[-X_a_Est1(1:3,iCount);X_a_Est1(4,iCount)];
   temp2 = fn_CrossTensor(X_a_true(1:4,iCount),0)*[-X_a_Est2(1:3,iCount);X_a_Est2(4,iCount)];
   err1(iCount) = asind(norm(temp1(1:3)));
   err2(iCount) = asind(norm(temp2(1:3)));
end
figure;
plot(err1);hold all;
plot(err2);
ylabel('$d\phi$ angle error(deg)','interpreter','latex','FontSize', 15);
legend('small-angle','Gibbs-vector');
xlabel('samples');
% for iCount = 1:4
%    subplot(2,2,iCount);   
%    plot(1:length(X_a_Est2),zeros(length(X_a_Est2),1),'k','linewidth',2,'linestyle','--');hold all;
%    %plot(X_a_Est1(iCount,:)-X_a_true(iCount,:),'linewidth',2);
%    %plot(X_a_Est2(iCount,:)-X_a_true(iCount,:),'linewidth',2);      
%    xlim([0,500]);
%    %ylim([-0.1,0.1]);
%    l = legend('zero line','small angle /delta','gibbs vector /delta');   
% end


% figure;
% title('omega and error details');
% for iCount = 1:3
%    subplot(3,2,2*iCount-1);
%    plot(X_a_true(iCount+4,:),'linewidth',2,'linestyle','-.');hold all;
%    plot(X_a_Est1(iCount+4,:),'linewidth',2,'linestyle','--');
%    plot(X_a_Est2(iCount+4,:),'linewidth',2,'linestyle',':');   
%    xlim([0,500]);   
%    legend('true','small angle','gibbs vector');   
% end
% for iCount = 1:3
%    subplot(3,2,2*iCount);   
%    plot(1:length(X_a_Est2),zeros(length(X_a_Est2),1),'k','linewidth',2,'linestyle','--');hold all;
%    plot(X_a_Est1(iCount+4,:)-X_a_true(iCount+4,:),'linewidth',2);
%    plot(X_a_Est2(iCount+4,:)-X_a_true(iCount+4,:),'linewidth',2);      
%    xlim([0,500]);
%    %ylim([-0.2,0.2]);
%    legend('zero line','small angle /delta','gibbs vector /delta');   
% end
% 
% figure;
% title('r and error details');
% for iCount = 1:3
%    subplot(3,2,2*iCount-1);
%    plot(X_a_true(iCount+7,:),'linewidth',2,'linestyle','-.');hold all;
%    plot(X_a_Est1(iCount+7,:),'linewidth',2,'linestyle','--');
%    plot(X_a_Est2(iCount+7,:),'linewidth',2,'linestyle',':');   
%    xlim([0,500]);   
%    legend('true','small angle','gibbs vector');   
% end
% for iCount = 1:3
%    subplot(3,2,2*iCount);   
%    plot(1:length(X_a_Est2),zeros(length(X_a_Est2),1),'k','linewidth',2,'linestyle','--');hold all;
%    plot(X_a_Est1(iCount+7,:)-X_a_true(iCount+7,:),'linewidth',2);
%    plot(X_a_Est2(iCount+7,:)-X_a_true(iCount+7,:),'linewidth',2);      
%    xlim([0,500]);
%    %ylim([-0.2,0.2]);
%    legend('zero line','small angle /delta','gibbs vector /delta');   
% end
% 
% figure;
% title('r_dot and error details');
% for iCount = 1:3
%    subplot(3,2,2*iCount-1);
%    plot(X_a_true(iCount+10,:),'linewidth',2,'linestyle','-.');hold all;
%    plot(X_a_Est1(iCount+10,:),'linewidth',2,'linestyle','--');
%    plot(X_a_Est2(iCount+10,:),'linewidth',2,'linestyle',':');   
%    xlim([0,500]);  
%    %ylim([-0.1,0.1]);
%    legend('true','small angle','gibbs vector');   
% end
% for iCount = 1:3
%    subplot(3,2,2*iCount);   
%    plot(1:length(X_a_Est2),zeros(length(X_a_Est2),1),'k','linewidth',2,'linestyle','--');hold all;
%    plot(X_a_Est1(iCount+10,:)-X_a_true(iCount+10,:),'linewidth',2);
%    plot(X_a_Est2(iCount+10,:)-X_a_true(iCount+10,:),'linewidth',2);      
%    xlim([0,500]);
%    ylim([-0.01,0.01]);
%    legend('zero line','small angle /delta','gibbs vector /delta');   
% end
% 
% 
% 
% 
% error_1 = X_a_true - X_a_Est1;
% error_2 = X_a_true - X_a_Est2;
% 
% sq_error_1 = error_1.*error_1;
% sq_error_2 = error_2.*error_2;
% 
% sum_sq_error_1 = sum(sq_error_1,1);
% sum_sq_error_2 = sum(sq_error_2,1);
% 
% figure;
% title('sum of squared errors');
% plot(1:length(X_a_Est2),-4*ones(length(X_a_Est2),1),'k','linewidth',2,'linestyle','--');hold all;
% plot(log(sum_sq_error_1)','linewidth',2);
% plot(log(sum_sq_error_2)','linewidth',2);
% xlim([0,500]);
% ylabel('log(sum of squared error)');
% legend('1e-3 line','small angle','gibbs vector');