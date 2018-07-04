function fig = checkVoltagePulse( sparams )
%CHECKVOLTAGEPULSE Summary of this function goes here
%   Detailed explanation goes here

    fig = figure;
    hold on;
    for ii = 1:sparams.numOfGates
        plot(sparams.voltagePulse(ii,:),'Linewidth',1.5);
    end
end

