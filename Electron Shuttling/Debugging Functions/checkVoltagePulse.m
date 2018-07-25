function fig = checkVoltagePulse( sparams )
%CHECKVOLTAGEPULSE Summary of this function goes here
%   Detailed explanation goes here

    fig = figure;
    hold on;
    set(gca,'Fontsize',14);
    xlabel('t/T','Interpreter','Latex','Fontsize',22);
    ylabel('Gate Voltage [V]','Interpreter','Latex','Fontsize',22);
    
    for ii = 1:sparams.numOfGates
        plot(linspace(0,1,length(sparams.voltagePulse(ii,:))),...
            sparams.voltagePulse(ii,:),'Linewidth',2.5);
    end
    ylim([min(min(sparams.voltagePulse(1,:))),max(max(sparams.voltagePulse(1,:)))*1.02])
    legend(sparams.gateLabels);
end

