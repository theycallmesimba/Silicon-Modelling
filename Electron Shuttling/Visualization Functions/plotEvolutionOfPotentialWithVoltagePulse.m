function fig = plotEvolutionOfPotentialWithVoltagePulse( sparams, xx )
%PLOTEVOLUTIONOFPOTENTIALWITHVOLTAGEPULSE Summary of this function goes here
%   Detailed explanation goes here

    % Now we need to make the individual gate interpolants for the pulse
    xPoints = linspace(0,1,length(sparams.voltagePulse(1,:)));
    qPoints = linspace(0,1,1001);
    vPulseGInterpolants = {};
    for vv = 1:sparams.numOfGates
        vPulseGInterpolants{vv} = griddedInterpolant({xPoints},sparams.voltagePulse(vv,:));
    end
    
    currPulse = getInterpolatedPulseValues(sparams,qPoints,vPulseGInterpolants);
    
    figure;
    hold on;
    
    xlabel('Position [nm]','Interpreter','Latex');
    for ii = 1:length(qPoints)
        currPot = sparams.P2DEGInterpolant(getInterpolantArgument(currPulse(:,ii),xx));
        currPot = squeezeFast(sparams.numOfGates,currPot)';
        if ii == 1
            plot(xx*1E-9,currPot);
            pause(0.01);
        end
        
        % Update the plot
        h = findobj('type','line');
        h(1).YData = currPot;
        pause(0.05);
    end
    delete(fig);
end

