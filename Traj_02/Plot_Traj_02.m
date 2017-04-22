figure;
for iCount = 1:4
    subplot(2,2,iCount);
    plot(t,eval(strcat('q_',num2str(iCount))),'-.','LineWidth',2); hold all;
    plot(t,eval(strcat('q_',num2str(iCount),'_ex')),'LineWidth',2);
    grid on;
    xlabel('time(sec)');    
    ylabel(strcat('$q_',num2str(iCount),'$'),'interpreter','latex','FontSize',20);
    legend('Simulated','Experimental');
end

