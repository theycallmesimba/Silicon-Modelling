function [vGTargMax, tcMax] = findResonantTunnelCoupling(sparams, xx, vVec, gateIndSweep, plotResult)
    % We want to find the point where the tunnel coupling is maximized
    % between two gates.  This is done by finding the voltage where the
    % gate being swept has maximum peak heights in both dots.
    
    % Need to find a good guess for where to look for the max tunnel
    % coupling
    % Just in case it's not.. Let's make our swept gate index be the same
    % voltage as the dot we will be coupled to.
    vVecpt1 = vVec;
    vVecpt1(gateIndSweep) = max(vVec);
    vVecCellpt1 = [num2cell(vVecpt1),xx];
    currPot = squeezeFast(sparams.numOfGates,sparams.P2DEGInterpolant(vVecCellpt1));
    
    [epsLpt, epsRpt] = getDetuning(sparams, xx, currPot.');
    
    % Now change the swept index by +/-.1 mV and see how detuning changes
    h = 0.1E-3;
    vVecpt2 = vVec;
    vVecpt2(gateIndSweep) = max(vVec) - h;
    vVecCellpt2 = [num2cell(vVecpt2),xx];
    currPot = squeezeFast(sparams.numOfGates,sparams.P2DEGInterpolant(vVecCellpt2));
    [epsLmh, epsRmh] = getDetuning(sparams, xx, currPot.');
    
    vVecpt3 = vVec;
    vVecpt3(gateIndSweep) = max(vVec) + h;
    vVecCellpt3 = [num2cell(vVecpt3),xx];
    currPot = squeezeFast(sparams.numOfGates,sparams.P2DEGInterpolant(vVecCellpt3));
    [epsLph, epsRph] = getDetuning(sparams, xx, currPot.');

    % Now find which voltage makes them 0 detuning
    temp1 = (epsLpt - epsLmh)/h;
    temp2 = (epsLph - epsLpt)/h;
    deL_dVi = (temp1 + temp2)/2;
    temp1 = (epsRpt - epsRmh)/h;
    temp2 = (epsRph - epsRpt)/h;
    deR_dVi = (temp1 + temp2)/2;
    
    eL0 = epsLpt - (deL_dVi*vVecpt1(gateIndSweep));
    eR0 = epsRpt - (deR_dVi*vVecpt1(gateIndSweep));
    
    vGuess = (eL0 - eR0)/(deR_dVi - deL_dVi);
    
    % Convert vVec from a list to a cell array and add on xx
    vVecCell = [num2cell(vVec),xx];
    
    deps = 5E-6*sparams.ee; % Have the searched voltage range be 10 ueV of detuning
    dVRange = abs(deps/deR_dVi);
    vRangePow = floor(log10(dVRange));
    options = optimset('TolX',10^(vRangePow - 2));
    vGTargMax = fminbnd(@(voltage) findPeakDifference(voltage),...
        vGuess - dVRange, vGuess + dVRange, options);
        
    function diff = findPeakDifference(currV)
        vVecCell{gateIndSweep} = currV;
        
        currPot = squeezeFast(sparams.numOfGates,sparams.P2DEGInterpolant(vVecCell));
        [currRho0, ~] = solve1DSingleElectronSE(sparams,1,xx,currPot);
        currRho0NSq = abs(currRho0).^2/norm(abs(currRho0).^2);
        
        pks = findpeaks(currRho0NSq);
        pks = pks(pks >= sparams.tcThreshold);
        
        if length(pks) == 2 
            diff = abs(pks(2) - pks(1));
        else
            diff = 1E3; % Some arbitrary large value
        end
    end
    
    vVecCell{gateIndSweep} = vGTargMax;
    currPotential = squeezeFast(sparams.numOfGates,sparams.P2DEGInterpolant(vVecCell));
    [currRhos, ens] = solve1DSingleElectronSE(sparams,2,xx,currPotential);

    currRho0NormSquared = abs(currRhos(:,1)').^2/norm(abs(currRhos(:,1)').^2);
    tcMax = abs(ens(2,2) - ens(1,1))/2;
    
    if plotResult
        tcFigure = figure('pos',[0 0 650 550]);
        movegui(tcFigure,'northeast');
        drawnow;               
        title(sprintf('Max Tc Point - Swept Gate %d',gateIndSweep));
        yyaxis left
        plot(xx,currPotential);
        yyaxis right
        plot(xx,currRho0NormSquared);

        pause(1.5);
        delete(tcFigure);
    end
end