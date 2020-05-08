function pulseGInterpolants = makePulseInterpolants( sparams, gridPoints, pulse )
%MAKEPULSEINTERPOLANTS Summary of this function goes here
%   Detailed explanation goes here

    pulseGInterpolants = {};
    for vv = 1:length(sparams.gatesUsedInPulse)
        pulseGInterpolants{vv} = griddedInterpolant({gridPoints},pulse(vv,:));
    end
end

