function shuttleSimulationAnalysis( sparams )
%SHUTTLESIMULATIONANALYSIS Summary of this function goes here

    fidFigure = figure('pos',[0 0 1100 900],'GraphicsSmoothing','on');
    movegui(fidFigure,'center');
    subplot(20,1,8:20);
    hold on;

    labelFontSize = 28;
    axesFontSize = 20;
    legFontSize = 22;
    
    set(gca,'Fontsize',axesFontSize);
    set(gcf,'color','w');
    for ii = [1,6,12,15,18,23]
        currFid = sparams.fidelity(ii,:);
        currFid(currFid == 0) = [];
        plot(linspace(0,1,length(currFid)),currFid,'Linewidth',2.5,...
            'DisplayName',sprintf('%0.3f ns',sparams.totalTime(ii)*1E9));
    end
    xlabel('$t/T$','Interpreter','Latex','Fontsize',labelFontSize);
    ylabel('Fidelity','Interpreter','Latex','Fontsize',labelFontSize);
    
    leg = legend('show','Location','northeastoutside');
    set(leg,'Fontsize',legFontSize);
    grid on;
    title(leg,'Total time $T$','Interpreter','Latex','Fontsize',legFontSize);
    posVec1 = get(gca,'Position');
    
%     [A,map] = rgb2ind(frame2im(getframe(fidFigure)),256);
%     imwrite(A,map, 'C:\Users\bbuonaco\Desktop\NetworkArchitectureFigures\fidelity.jpg');

    subplot(20,1,1:6);
    hold on;
    posVec2 = get(gca,'Position');
    set(gca,'Fontsize',axesFontSize,'Position',[posVec2(1:2) posVec1(3)*0.867 posVec2(4)]);
    for ii = 1:sparams.numOfGates
        plot(linspace(0,1,length(sparams.voltagePulse(ii,:))),sparams.voltagePulse(ii,:),'Linewidth',2.5);
    end
    ylim([0.6,0.825]);
%     xlabel('Time [ns]','Interpreter','Latex','Fontsize',labelFontSize);
    ylabel('$V_g$ [V]','Interpreter','Latex','Fontsize',labelFontSize);
    leg = legend('V_1','V_2','V_3','Location','northeastoutside');
    set(leg,'Fontsize',legFontSize);
    set(gca,'xtick',[]);
end

