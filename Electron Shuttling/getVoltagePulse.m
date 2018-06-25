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
    vMinBnd = 0.786;
    vMaxBnd = 0.796;
    
    [g2TurnTcOn, g2Max] = findTunnelCouplingTurnOn(sparams, xx,...
        [g1Max,g2Min,g3Min,g4Min,g5Min], vMinBnd, vMaxBnd, 2);
    [g1TurnTcOff, ~] = findTunnelCouplingTurnOn(sparams, xx,...
        [g1Min,g2Max,g3Min,g4Min,g5Min], vMinBnd, vMaxBnd, 1);
    
    [g3TurnTcOn, g3Max] = findTunnelCouplingTurnOn(sparams, xx,...
        [g1Min,g2Max,g3Min,g4Min,g5Min], vMinBnd, vMaxBnd, 3);
    [g2TurnTcOff, ~] = findTunnelCouplingTurnOn(sparams, xx,...
        [g1Min,g2Min,g3Max,g4Min,g5Min], vMinBnd, vMaxBnd, 2);
    
    [g4TurnTcOn, g4Max] = findTunnelCouplingTurnOn(sparams, xx,...
        [g1Min,g2Min,g3Max,g4Min,g5Min], vMinBnd, vMaxBnd, 4);
    [g3TurnTcOff, ~] = findTunnelCouplingTurnOn(sparams, xx,...
        [g1Min,g2Min,g3Min,g4Max,g5Min], vMinBnd, vMaxBnd, 3);
    
    [g5TurnTcOn, g5Max] = findTunnelCouplingTurnOn(sparams, xx,...
        [g1Min,g2Min,g3Min,g4Max,g5Min], vMinBnd, vMaxBnd, 5);
    [g4TurnTcOff, ~] = findTunnelCouplingTurnOn(sparams, xx,...
        [g1Min,g2Min,g3Min,g4Min,g5Max], vMinBnd, vMaxBnd, 4);
    g5Max = 0.8;
    
%     g2Max = findZeroDetuning(sparams,xx,[g1Max,g2Min,g3Min,g4Min,g5Min],2);
%     g2Max = 0.7929;
%     g3Max = 0.7927;
%     g3Max = findZeroDetuning(sparams,xx,[g1Min,g2Max,g3Min,g4Min,g5Min],3);
%     g4Max = 0.7929;
%     g4Max = findZeroDetuning(sparams,xx,[g1Min,g2Min,g3Min,g4Min,g5Max],4);
%     ratio = 0.989;
        
    gPulse = {};
    
    gPulse{1,1} = [g1Max, g1Max, g1TurnTcOff, g1Min, g1Min];
    gPulse{1,2} = [0, 12.5, 24, 25, 100];
    
    gPulse{2,1} = [g2Min, g2TurnTcOn, g2Max, g2Max, g2TurnTcOff, g2Min, g2Min];
    gPulse{2,2} = [0, 1, 12.5, 37.5, 49, 50, 100];
    
    gPulse{3,1} = [g3Min, g3Min, g3TurnTcOn, g3Max, g3Max, g3TurnTcOff, g3Min, g3Min];
    gPulse{3,2} = [0, 25, 26, 37.5, 62.5, 74, 75, 100];
    
    gPulse{4,1} = [g4Min, g4Min, g4TurnTcOn, g4Max, g4Max, g4TurnTcOff, g4Min];
    gPulse{4,2} = [0, 50, 51, 62.5, 87.5, 99, 100];
    
    gPulse{5,1} = [g5Min, g5Min, g5TurnTcOn, g5Max, g5Max];
    gPulse{5,2} = [0, 75, 76, 87.5, 100];
    
    for ii = 1:sparams.numOfGates
        sparams.voltagePulse(ii,:) = interp1(gPulse{ii,2},gPulse{ii,1},0:100);
    end
end

function [vGTargTurnOn, vGTargMax] = findTunnelCouplingTurnOn(sparams, xx, vVec,...
    vMinBnd, vMaxBnd, gateIndSweep)
    % What we want to do here is find the point in time where tunnel
    % coupling turns on.  So we will slowly sweep through the voltage (in
    % 0.01 mV increments until we find 2 peaks of the wavefunction (one of
    % the peaks should be incredibly small).  Then we will continue
    % sweeping until the peak heights are close to equal indicating maximal
    % tunnel coupling.  The first voltage point will be where we want to
    % sweep fast to and from there on we will sweep slowly.

    % Convert vVec from a list to a cell array and add on xx
    vVecCell = [num2cell(vVec),xx];
    
    turnOnFound = 0;
    maxTcFound = 0;
    
    currPeakDifference = 1E10;
    tcFigure = figure('pos',[0 0 1300 550]);
    movegui(tcFigure,'northeast');
    drawnow;
    title(sprintf('Finding Tc Points for Gate = %d',gateIndSweep));
    
    for currV = vMinBnd:1E-5:vMaxBnd
        if maxTcFound
            break;
        end
        vVecCell{gateIndSweep} = currV;
        currPotential = squeeze(sparams.P2DEGInterpolant(vVecCell));
        [currRho0, ~] = solve1DSingleElectronSE(sparams,1,xx,currPotential);
        currRho0NormSquared = abs(currRho0).^2/norm(abs(currRho0).^2);
        
        pks = findpeaks(currRho0NormSquared);
        pks = pks(pks >= 0.0001);
                
        if length(pks) == 2 
            if ~turnOnFound
                turnOnFound = 1;
                vGTargTurnOn = currV;

                findpeaks(currRho0NormSquared);

                subplot(1,2,1);
                title(sprintf('Tc Turn On - Gate %d',gateIndSweep));
                yyaxis left
                plot(xx,currPotential);
                yyaxis right
                plot(xx,currRho0NormSquared);
                drawnow;
            end
            
            diff = abs(pks(2) - pks(1));
            if currPeakDifference > diff
                currPeakDifference = diff;
                vGTargMax = currV;
            else
                maxTcFound = 1;
                
                subplot(1,2,2);
                title(sprintf('Max Tc - Gate %d',gateIndSweep));
                yyaxis left
                plot(xx,currPotential);
                yyaxis right
                plot(xx,currRho0NormSquared);
            end
        end
    end
    pause(4);
    delete(tcFigure);
end

function vGTarg = findZeroDetuning( sparams, xx, vVec, gateIndSweep )
%Finds what gate voltage value will give you zero detuning
%between two gates.  The voltage on one gate must be fixed and the other
%gate voltage will vary until 0 detuning is achieved.  It is assumed you
%are choosing a fixed gate based on the inputted voltage vector (vVec)

    % Convert vVec from a list to a cell array and add on xx
    vVecCell = [num2cell(vVec),xx];
    
    % Get voltage range for that gate
    gateValues = reshape([sparams.potentials.gateValues],sparams.numOfGates,[])';
    minGVal = min(unique(gateValues(:,gateIndSweep)));
    maxGVal = max(unique(gateValues(:,gateIndSweep)));
    
    % Perform the minimization
    scaleFactor = 1000;
    options = optimset('TolX',1E-6,'Display','iter');
    [vGTarg, ~] = fminbnd(@(x) findMinDeltaE(x),minGVal*scaleFactor,maxGVal*scaleFactor,options);
    vGTarg = vGTarg/scaleFactor;
    function deltaE = findMinDeltaE(gateValue)
        % Update the current voltage vector and get the corresponding
        % potential
        vVecCell{gateIndSweep} = gateValue/scaleFactor;
        currPot = squeeze(sparams.P2DEGInterpolant(vVecCell));
        
        % Find the wells of the potential calculate their detuning
        peaks = sort(findpeaks(-currPot),'descend');
        deltaE = (peaks(1) - peaks(2))/sparams.ee;
    end
end







