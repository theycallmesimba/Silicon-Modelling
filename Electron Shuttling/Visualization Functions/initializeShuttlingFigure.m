function fig = initializeShuttlingFigure( sparams, vPulse, wf1, wf2, xx, time )
%INITIALIZESHUTTLINGFIGURE Summary of this function goes here
%   Detailed explanation goes here


    vvInitial = sparams.P2DEGInterpolant([num2cell(vPulse(:,1)'),xx]);
    vvInitial = squeezeFast(length(sparams.gatesUsedInPulse),vvInitial);
    
    fig = figure('pos',[0 0 1300 650],'Color','white');
    % Set up colors
    cm = viridis();
    set(fig,'defaultAxesColorOrder',[cm(20,:); cm(110,:);]);
    movegui(fig,'northeast');
    
    xx = xx/1E-9;
    axisFontSize = 20;
    tickLabelFontSize = 20;
    labelFontSize = 35;
    titleFontSize = 35;
    legenedFontSize = 25;
    
    hold on;
    set(gca,'TickLabelInterpreter','latex','Fontsize',tickLabelFontSize);
    set(gca,'Fontsize',axisFontSize);
%     title(['Shuttling Simulation ' num2str(time) '[s]'],'interpreter','latex','fontsize',titleFontSize);
    xlabel('Position [nm]','Interpreter','Latex','Fontsize',labelFontSize);
    xlim([min(xx),max(xx)]);
    
    yyaxis left
%     CBoffset = -2.2747;
%     vvTemp = vvInitial/sparams.ee - CBoffset;
    plot(xx,vvInitial/sparams.ee,'Linewidth',2);
    plot(xx,vvTemp,'Linewidth',2);
    ylabel('Energy [eV]','Interpreter','Latex','Fontsize',labelFontSize);
    
    yyaxis right
    plot(xx,abs(wf1).^2/norm(abs(wf1).^2),'Linewidth',2);
    plot(xx,abs(wf2).^2/norm(abs(wf2).^2),'Linewidth',2);
    ylabel('Probability','Interpreter','Latex','Fontsize',labelFontSize);
    
    legend({'${\rm Potential}$','$|\psi_{\rm Sim}|^2$','$|\psi_0|^2$'},...
        'Interpreter','latex','Fontsize',legenedFontSize);
end

