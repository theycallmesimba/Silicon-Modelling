function sparams = makePotentialsInterpolantsAlt( sparams, xx, yy )
%MAKEPOTENTIALSINTERPOLANTS Summary of this function goes here
%   Detailed explanation goes here

    % Get all the gate values loaded
    allGateValues = reshape([sparams.potentials.gateValues],sparams.numOfGates,[])';

    % Make the function argument for the grid vectors for loadable gate
    % voltages which is used to make the gridded interpolant
    argGridVecVolts = {};
    for ii = 1:length(sparams.interpableGates)
        argGridVecVolts = [argGridVecVolts,unique(allGateValues(:,sparams.interpableGates(ii)))];
    end

    % Make the function argument for the magnitude of the grid vectors
    % which is used to reshape the structure variable
    % sparams.potentials.pot2D(EG)
    argReshapeMagnitudes = [sparams.nyGrid,sparams.nxGrid];
    for ii = length(sparams.interpableGates):-1:1
        argReshapeMagnitudes = [argReshapeMagnitudes,length(argGridVecVolts{ii})];
    end

    % Need to permute the reshaped sparams.potentials.pot2D(EG)
    argPermutation2D = [(length(sparams.gatesUsedInPulse):-1:1) + 2,1,2];
%     argPermutation2DEG = [(length(sparams.gatesUsedInPulse):-1:1) + 1,1];
    
    % Create the potentials interpolant
    sparams.P2DInterpolant = griddedInterpolant([argGridVecVolts,yy,xx],...
        permute(reshape([sparams.potentials.pot2D],...
        argReshapeMagnitudes),argPermutation2D),sparams.interpType,sparams.extrapType);

%     sparams.P2DEGInterpolant = griddedInterpolant([argGridVecVolts,xx],...
%         permute(reshape([sparams.potentials.pot2DEG],...
%         argReshapeMagnitudes(2:end)),argPermutation2DEG),sparams.interpType,sparams.extrapType);
end

