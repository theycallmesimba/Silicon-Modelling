function sparams = getVoltagePulse( sparams, xx )
%GETVOLTAGEPULSE Summary of this function goes here
%   Detailed explanation goes here

    % The idea behind this function is to generate a pulsing sequence for
    % all the gates we have control over in our geometry.
    % We simply define each pulsing sequence according to percentage of the
    % total time.  For now, we use a grid of 101 (+1 for t0) points.  So we have
    % accuracy of our control pulses to within 1% of our total time.  This
    % is easily adjustable if we need in the future.
    
    sparams.voltagePulse = zeros(sparams.numOfGates,101);
    
    g1Min = 0.6;
    g2Min = 0.6;
    g3Min = 0.6;
    g4Min = 0.6;
    g5Min = 0.6;
    
    g1Max = 0.8;
    g5Max = 0.8;
%     g2Max = findZeroDetuning(sparams,xx,[g1Max,g2Min,g3Min,g4Min,g5Min],2);
    g2Max = 0.7929;
%     plotPotentialAndGroundWF(sparams,[g1Max,g2Max,g3Min,g4Min,g5Min],xx);
    g3Max = 0.7927;
%     g3Max = findZeroDetuning(sparams,xx,[g1Min,g2Max,g3Min,g4Min,g5Min],3);
%     plotPotentialAndGroundWF(sparams,[g1Min,g2Max,g3Max,g4Min,g5Min],xx);
    g4Max = 0.7929;
%     g4Max = findZeroDetuning(sparams,xx,[g1Min,g2Min,g3Min,g4Min,g5Max],4);
%     plotPotentialAndGroundWF(sparams,[g1Min,g2Min,g3Min,g4Max,g5Max],xx);
    
    ratio = 0.989;
        
    gPulse = {};
    
    gPulse{1,1} = [g1Max, g1Max, ratio*g1Max, g1Min, g1Min];
    gPulse{1,2} = [0, 12.5, 24, 25, 100];
    
    gPulse{2,1} = [g2Min, ratio*g2Max, g2Max, g2Max, ratio*g2Max, g2Min, g2Min];
    gPulse{2,2} = [0, 1, 12.5, 37.5, 49, 50, 100];
    
    gPulse{3,1} = [g3Min, g3Min, ratio*g3Max, g3Max, g3Max, ratio*g3Max, g3Min, g3Min];
    gPulse{3,2} = [0, 25, 26, 37.5, 62.5, 74, 75, 100];
    
    gPulse{4,1} = [g4Min, g4Min, ratio*g4Max, g4Max, g4Max, ratio*g4Max, g4Min];
    gPulse{4,2} = [0, 50, 51, 62.5, 87.5, 99, 100];
    
    gPulse{5,1} = [g5Min, g5Min, ratio*g5Max, g5Max, g5Max];
    gPulse{5,2} = [0, 75, 76, 87.5, 100];
    
    for ii = 1:sparams.numOfGates
        sparams.voltagePulse(ii,:) = interp1(gPulse{ii,2},gPulse{ii,1},0:100);
    end
end











