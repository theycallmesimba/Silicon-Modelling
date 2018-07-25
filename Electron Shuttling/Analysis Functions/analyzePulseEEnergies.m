function analyzePulseEEnergies( sparams, xx )
%ANALYZEPULSEEIEGENENERGIES Summary of this function goes here
%   Detailed explanation goes here

    sparams.vPulseGInterpolants = {};
    xPoints = linspace(0,1,length(sparams.voltagePulse(1,:)));    

    for vv = 1:sparams.numOfGates
        vPulseGInterpolants{vv} = griddedInterpolant({xPoints},sparams.voltagePulse(vv,:));
    end
    
    qPots = linspace(0,1,301);
    pulse = getInterpolatedPulseValues(sparams,qPots,vPulseGInterpolants);
    
    nSols = 4;
    ens = zeros(length(qPots),nSols);
    for ii = 1:length(qPots)
        currPot = sparams.P2DEGInterpolant(getInterpolantArgument(pulse(:,ii),xx));
        [~,energies] = solve1DSingleElectronSE(sparams, nSols, xx, currPot);
        
        ens(ii,:) = diag(energies)'/sparams.ee;
    end
    
    plotFunctionOverVoltagePulse(sparams,xPoints,qPots,ens);
end

