function sparams = makePotentialsInterpolants( sparams, xx, zz )
%MAKEPOTENTIALSINTERPOLANTS Summary of this function goes here
%   Detailed explanation goes here

    % Get the min and max value for each gate
    allGateValues = reshape([sparams.potentials.gateValues],sparams.numOfGates,[])';

    argGridVecVolts = {};
    for ii = 1:sparams.numOfGates
        argGridVecVolts = [argGridVecVolts,unique(allGateValues(:,ii))];
    end
    
    argReshapeMagnitudes = [sparams.nzGrid,sparams.nxGrid];
    for ii = sparams.numOfGates:-1:1
        argReshapeMagnitudes = [argReshapeMagnitudes,length(argGridVecVolts{ii})];
    end
    
    argPermutation2D = [(sparams.numOfGates:-1:1) + 2,1,2];
    argPermutation2DEG = [(sparams.numOfGates:-1:1) + 1,1];
    
    % Create the potentials interpolant
    sparams.P2DInterpolant = griddedInterpolant([argGridVecVolts,zz,xx],...
        permute(reshape([sparams.potentials.pot2D],...
        argReshapeMagnitudes),argPermutation2D));
 
    sparams.P2DEGInterpolant = griddedInterpolant([argGridVecVolts,xx],...
        permute(reshape([sparams.potentials.pot2DEG],...
        argReshapeMagnitudes(2:end)),argPermutation2DEG));
end

