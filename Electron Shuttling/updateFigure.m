function updateFigure( sparams, fig, wf1, wf2, xx, vv, timeInd )
%UPDATEFIGURE Summary of this function goes here
%   Detailed explanation goes here

    figure(fig);
    clf;
    hold on;
    plot(xx,vv/sparams.ee);
    plot(xx,abs(wf1.^2)/2000000000 + sparams.figWFMin/sparams.ee);
    plot(xx,abs(wf2.^2)/2000000000 + sparams.figWFMin/sparams.ee);
    title(['Shuttling Simulation ' num2str(sparams.totalTime(timeInd)) '[s]'],'interpreter','latex','fontsize',12);
    xlabel('Position [m]','interpreter','latex','fontsize',12);
    ylabel('Energy [eV]','interpreter','latex','fontsize',12);
    legend('Current Potential','Current Sim \psi','Current Ground State \psi');
end

