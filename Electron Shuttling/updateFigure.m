function updateFigure( sparams, fig, wf1, wf2, xx, vv, timeInd )
%UPDATEFIGURE Summary of this function goes here
%   Detailed explanation goes here

    figure(fig);
    clf;
    hold on;
    title(['Shuttling Simulation ' num2str(sparams.totalTime(timeInd)) '[s]'],'interpreter','latex','fontsize',16);
    xlabel('Position [nm]','Interpreter','Latex','Fontsize',14);
    xlim([min(xx),max(xx)]);
    
    yyaxis left
    plot(xx,vv/sparams.ee,'Linewidth',1.5);
    ylabel('Potential [V]','Interpreter','Latex','Fontsize',14);
    
    yyaxis right
    plot(xx,abs(wf1).^2/norm(abs(wf1).^2),'Linewidth',1.5);
    plot(xx,abs(wf2).^2/norm(abs(wf2).^2),'Linewidth',1.5);
    ylabel('Probability [Arb]','Interpreter','Latex','Fontsize',14);
    
    legend('Current Potential','Current Sim |\psi|^2','Current Ground State |\psi|^2');
end

