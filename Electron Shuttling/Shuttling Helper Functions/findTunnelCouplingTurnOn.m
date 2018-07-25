function vGTargTurnOn = findTunnelCouplingTurnOn(sparams, xx, vVec, vMinBnd,...
    vMaxBnd, gateIndSweep)
    % We want to find the point where the tunnel coupling is just turn on
    % between two gates.  We do this by doing a minimization over a
    % function that calculates the peak heights of the wave function for a
    % given voltage and seeing when there just begins to be a peak in the
    % swept gate index potential well.

    % Convert vVec from a list to a cell array and add on xx
    vVecCell = [num2cell(vVec),xx];
    
    options = optimset('TolX',1E-5);
    vGTargTurnOn = fminbnd(@(voltage) findPeakHeight(voltage),...
        vMinBnd, vMaxBnd, options);
        
    function diffThresh = findPeakHeight(currV)
        vVecCell{gateIndSweep} = currV;
        
        currPot = squeezeFast(sparams.numOfGates,sparams.P2DEGInterpolant(vVecCell));
        [currRho0, ~] = solve1DSingleElectronSE(sparams,1,xx,currPot);
        currRho0NSq = abs(currRho0).^2/norm(abs(currRho0).^2);
        
        pks = findpeaks(currRho0NSq);
        pks = pks(pks >= sparams.tcTuningThreshold);
        
        if length(pks) == 2
            pks = sort(pks);
            diffThresh = pks(1) - sparams.tcTuningThreshold;
        else
            diffThresh = 1E3;
        end
    end
    
    vVecCell{gateIndSweep} = vGTargTurnOn;
    currPotential = squeezeFast(sparams.numOfGates,sparams.P2DEGInterpolant(vVecCell));
    [currRhos, ~] = solve1DSingleElectronSE(sparams,2,xx,currPotential);

    currRho0NormSquared = abs(currRhos(:,1)').^2/norm(abs(currRhos(:,1)').^2);
    
    tcFigure = figure('pos',[0 0 650 550]);
    movegui(tcFigure,'northeast');
    drawnow;               
    title(sprintf('Tc Turn On/Off Point - Swept Gate %d',gateIndSweep));
    yyaxis left
    plot(xx,currPotential);
    yyaxis right
    plot(xx,currRho0NormSquared);
    
    pause(2);
    delete(tcFigure);
end