load('data_trajectory_T6_45\Data_01.mat');
for iCount = 1:4
   subplot(2,2,iCount);
   plot(Mu(:,iCount));
end