function detInterpolants = makeDetuningInterpolants(gridPoints, epsL, epsR)
%GETDETUNINGINTERPOLANTS Summary of this function goes here
%   Detailed explanation goes here
    
    detInterpolants = {};
    detInterpolants{1} = griddedInterpolant({gridPoints},epsL);
    detInterpolants{2} = griddedInterpolant({gridPoints},epsR);
end

