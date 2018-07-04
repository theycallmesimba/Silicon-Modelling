function plotFunctionOverVoltagePulse( sparams, xPulse, xxDataToPlot, yyDataToPlot )
%PLOTFUNCOVERVOLTAGEPULSE Summary of this function goes here
%   Detailed explanation goes here

    figure;
    axisFontSize = 20;
    tickFontSize = 15;
    titleFontSize = 25;
    set(gca,'Fontsize',tickFontSize);
    hold on;
%     title('E-Energies vs. Pulse','Interpreter','Latex','Fontsize',titleFontSize);
%     xlabel('\% of shuttle time','Interpreter','Latex','Fontsize',axisFontSize);
    
    yyaxis left
    plot(xxDataToPlot,yyDataToPlot,'Linewidth',2.0,'Linestyle','-');
%     ylabel('Energy [eV]','Interpreter','Latex','Fontsize',axisFontSize);
    
    yyaxis right
    for ii = 1:sparams.numOfGates
        plt = plot(xPulse,sparams.voltagePulse(ii,:),'Linewidth',2.0,...
            'Linestyle','-','Marker','none');
        plt.Color(4) = 0.7;
    end
    ylabel('Voltage [V]','Interpreter','Latex','Fontsize',axisFontSize);
end

