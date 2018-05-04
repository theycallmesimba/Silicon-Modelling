function ind = getPotentialGivenGateValues( sparams, gateVals )
%GETPOTENTIALGIVENGATEVALUES Summary of this function goes here
%   Detailed explanation goes here

    allGateValues = reshape([sparams.potentials.gateValues],sparams.numOfGates,[])';
    [~,ind] = ismembertol(gateVals,allGateValues,'ByRows',true);
end

