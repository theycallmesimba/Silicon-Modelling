function plot2DBasis( gridparams, basisToPlot, nStates )
%PLOT2DBASIS Summary of this function goes here
%   Detailed explanation goes here

    XX = gridparams.XX;
    YY = gridparams.YY;

    maxPlotsPerFigure = 8;

    for ii = 1:nStates
        currWFMG = basisToPlot(ii).wavefunctionMG;
       
        if mod(ii-1,maxPlotsPerFigure) == 0
            fig = figure('Color','white');
            pos = get(fig,'position');
            set(fig,'position',[pos(1:2)/4 pos(3)*1.8 pos(4)*1.3]);
        end
        subplot(2,maxPlotsPerFigure/2,mod(ii-1,maxPlotsPerFigure)+1);
        
        s = surf(XX,YY,currWFMG);
        set(s,'edgecolor','none');
        set(gca,'Fontsize',10,'TicklabelInterpreter','latex');
        xlabel('x-axis','Fontsize',14,'Interpreter','latex');
        ylabel('y-axis','Fontsize',14,'Interpreter','latex');
        xlim([min(min(XX)),max(max(XX))]);
        ylim([min(min(YY)),max(max(YY))]);
        view(2);
        colormap(viridis());
        
        drawnow;
    end
end

