function updateFigure( sparams, fig, wf1, wf2, xx, vv, timeInd )
%UPDATEFIGURE Summary of this function goes here
%   Detailed explanation goes here

    xx = xx/1E-9;
    axisFontSize = 16;
    labelFontSize = 30;
    titleFontSize = 35;
    
    figure(fig);
    clf;
    hold on;
    set(gca,'Fontsize',axisFontSize);
    title(['Shuttling Simulation ' num2str(sparams.totalTime(timeInd)) '[s]'],'interpreter','latex','fontsize',titleFontSize);
    xlabel('Position [nm]','Interpreter','Latex','Fontsize',labelFontSize);
    xlim([min(xx),max(xx)]);
    
    yyaxis left
    plot(xx,vv/sparams.ee,'Linewidth',1.5);
    ylabel('Potential [V]','Interpreter','Latex','Fontsize',labelFontSize);
    
    yyaxis right
    plot(xx,abs(wf1).^2/norm(abs(wf1).^2),'Linewidth',2.5);
    plot(xx,abs(wf2).^2/norm(abs(wf2).^2),'Linewidth',2.5);
    ylabel('Probability','Interpreter','Latex','Fontsize',labelFontSize);
    
    legend('Potential','Sim |\psi|^2','Ground State |\psi|^2');
end

