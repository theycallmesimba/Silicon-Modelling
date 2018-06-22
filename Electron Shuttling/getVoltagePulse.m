function sparams = getVoltagePulse( sparams, xx )
%GETVOLTAGEPULSE Summary of this function goes here
%   Detailed explanation goes here

    % The idea behind this function is to generate a pulsing sequence for
    % all the gates we have control over in our geometry.
    % We simply define each pulsing sequence according to percentage of the
    % total time.  For now, we use a grid of 100 points.  So we have
    % accuracy of our control pulses to within 1% of our total time.  This
    % is easily adjustable if we need in the future.
    
    sparams.voltagePulse = zeros(sparams.numOfGates,101);

    % Let's do the first part of the sweep
    g1Min = 0.6;
    g1Max = 0.8;
    g2Min = 0.6;
    g3Min = 0.6;
    g3Max = 0.8;
    ratio = 0.989;
    sparams.gateLabels = {'V_{L1}','V_{L2}','V_{C}','V_{R2}','V_{R1}'};

%     g1Min = 1.5;
%     g1Max = 1.7;
%     g2Min = 1.5; 
%     g3Min = 1.5;
%     g3Max = 1.7;
%     ratio = 0.99;
    
    % Now, we wish to find what value the second gate needs to be so that
    % the tunnel coupling is maximal.  Due to cross capacitances, setting
    % V1 = V2 does not actually mean the detuning is 0.  This short segment
    % finds that value
    scaleFactor = 1000;
    [g2Max, ~] = fminbnd(@(x) findMinDeltaE(x),0.7*scaleFactor,0.81*scaleFactor);
    g2Max = g2Max/scaleFactor;
    function deltaE = findMinDeltaE(g2)
        g2 = g2/scaleFactor;
        currPot = squeeze(sparams.P2DEGInterpolant({g1Max,g2,g3Min,xx}));
            
        peaks = sort(findpeaks(-currPot),'descend');
        deltaE = (peaks(1) - peaks(2))/sparams.ee;         
    end
    
    ind = 1:3;
    gate1p = ones(1,length(ind))*g1Max;
    gate2p = linspace(g2Min,ratio*g2Max,length(ind));
    gate3p = ones(1,length(ind))*g3Min;
    
    ind = 3:26;
    temp1 = ones(1,length(ind))*g1Max;
    temp2 = linspace(ratio*g2Max,g2Max,length(ind));
    temp3 = ones(1,length(ind))*g3Min;
    gate1p = [gate1p, temp1(2:end)];
    gate2p = [gate2p, temp2(2:end)];
    gate3p = [gate3p, temp3(2:end)];
    
    ind = 26:49;
    temp1 = linspace(g1Max,ratio*g1Max,length(ind));
    temp2 = ones(1,length(ind))*g2Max;
    temp3 = ones(1,length(ind))*g3Min;
    gate1p = [gate1p, temp1(2:end)];
    gate2p = [gate2p, temp2(2:end)];
    gate3p = [gate3p, temp3(2:end)];
    
    ind = 49:51;
    temp1 = linspace(ratio*g1Max,g1Min,length(ind));
    temp2 = ones(1,length(ind))*g2Max;
    temp3 = ones(1,length(ind))*g3Min;
    gate1p = [gate1p, temp1(2:end)];
    gate2p = [gate2p, temp2(2:end)];
    gate3p = [gate3p, temp3(2:end)];
    
    ind = 51:53;
    temp1 = ones(1,length(ind))*g1Min;
    temp2 = ones(1,length(ind))*g2Max;
    temp3 = linspace(g3Min,ratio*g3Max,length(ind));
    gate1p = [gate1p, temp1(2:end)];
    gate2p = [gate2p, temp2(2:end)];
    gate3p = [gate3p, temp3(2:end)];
    
    ind = 53:76;
    temp1 = ones(1,length(ind))*g1Min;
    temp2 = ones(1,length(ind))*g2Max;
    temp3 = linspace(ratio*g3Max,g3Max,length(ind));
    gate1p = [gate1p, temp1(2:end)];
    gate2p = [gate2p, temp2(2:end)];
    gate3p = [gate3p, temp3(2:end)];
    
    ind = 76:98;
    temp1 = ones(1,length(ind))*g1Min;
    temp2 = linspace(g2Max,ratio*g2Max,length(ind));
    temp3 = ones(1,length(ind))*g3Max;
    gate1p = [gate1p, temp1(2:end)];
    gate2p = [gate2p, temp2(2:end)];
    gate3p = [gate3p, temp3(2:end)];
    
    ind = 98:100;
    temp1 = ones(1,length(ind))*g1Min;
    temp2 = linspace(ratio*g2Max,g2Min,length(ind));
    temp3 = ones(1,length(ind))*g3Max;
    gate1p = [gate1p, temp1(2:end)];
    gate2p = [gate2p, temp2(2:end)];
    gate3p = [gate3p, temp3(2:end)];
    
    gate1p = [gate1p, g1Min];
    gate2p = [gate2p, g2Min];
    gate3p = [gate3p, g3Max];
    
    sparams.voltagePulse(1,:) = gate1p;
    sparams.voltagePulse(2,:) = gate2p;
    sparams.voltagePulse(3,:) = gate3p;
    
end











