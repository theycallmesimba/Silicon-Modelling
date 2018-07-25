function [vGTargMax, tcMax] = findTunnelCouplingMax(sparams, xx, vVec, vMinBnd,...
    vMaxBnd, gateIndSweep)
    % We want to find the point where the tunnel coupling is maximized
    % between two gates.  This is done by finding the voltage where the
    % gate being swept has maximum peak heights in both dots.

    % Convert vVec from a list to a cell array and add on xx
    vVecCell = [num2cell(vVec),xx];
    
    options = optimset('TolX',1E-5);
    vGTargMax = fminbnd(@(voltage) findPeakDifference(voltage),...
        vMinBnd, vMaxBnd, options);
        
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
    tcMax = abs(ens(2,2) - ens(1,1));
    
    tcFigure = figure('pos',[0 0 650 550]);
    movegui(tcFigure,'northeast');
    drawnow;               
    title(sprintf('Max Tc Point - Swept Gate %d',gateIndSweep));
    yyaxis left
    plot(xx,currPotential);
    yyaxis right
    plot(xx,currRho0NormSquared);
    
    pause(1);
    delete(tcFigure);
end