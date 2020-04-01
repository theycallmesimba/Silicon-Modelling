function plotFunctionOverVoltagePulse( sparams, xxData, yDataPulse, yyDataToPlot )
%PLOTFUNCOVERVOLTAGEPULSE Summary of this function goes here
%   Detailed explanation goes here

    figure;
    axisFontSize = 25;
    tickFontSize = 18;

    set(gca,'Fontsize',tickFontSize);
    hold on;
%     title('E-Energies vs. Pulse','Interpreter','Latex','Fontsize',titleFontSize);
%     xlabel('\% of shuttle time','Interpreter','Latex','Fontsize',axisFontSize);
    
    yyaxis right
    for ii = 1:sparams.numOfGates
        plt = plot(xxData,yDataPulse(ii,:),'Linewidth',2.0,...
            'Linestyle','-','Marker','none');
    end
    ylabel('Voltage [V]','Interpreter','Latex','Fontsize',axisFontSize);

    yyaxis left
    plot(xxData,yyDataToPlot,'Linewidth',2.0,'Linestyle','-');
    xlim([min(xxData),max(xxData)]);
    xlabel('Time [s]','Interpreter','Latex');
    ylabel('Adiabatic Parameter','Interpreter','Latex');
   
end

