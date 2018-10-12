function detInterpolants = makeDetuningInterpolants(gridPoints, effHamiltonianParams)
%GETDETUNINGINTERPOLANTS Summary of this function goes here
%   Detailed explanation goes here

    [epsL, epsR, ~, ~, ~, ~, ~, ~] = decodeEffHamiltonianParamVariable(effHamiltonianParams);
    
    detInterpolants = {};
    detInterpolants{1} = griddedInterpolant({gridPoints},epsL);
    detInterpolants{2} = griddedInterpolant({gridPoints},epsR);
end

