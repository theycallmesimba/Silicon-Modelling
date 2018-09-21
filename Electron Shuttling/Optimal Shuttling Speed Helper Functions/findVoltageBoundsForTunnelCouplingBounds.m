function voltageBounds = findVoltageBoundsForTunnelCouplingBounds( sparams, xx, tcBounds, vSearchBnds, dotLocs)
%FINDVOLTAGEBOUNDSFORTUNNELCOUPLINGBOUNDS Summary of this function goes here
%   Detailed explanation goes here

    fprintf(1,'Finding voltages for desired tc boundary...\n');

    % 1 uV accuracy in finding the voltage bounds
    options = optimset('TolX',1E-5);
    
    % First get lower voltage bound (i.e. highest tc bound point)
    voltageBounds(1) = fminbnd(@(voltage) findBoundaryPoint(voltage, max(tcBounds)),...
         min(vSearchBnds), max(vSearchBnds), options);

    % Then get higher voltage bound (i.e. lowest tc bound point)
    voltageBounds(2) = fminbnd(@(voltage) findBoundaryPoint(voltage, min(tcBounds)),...
        min(vSearchBnds), max(vSearchBnds), options);
    
    function diff = findBoundaryPoint(currV, tcBound)
        [~, currTc] = findResonantTunnelCoupling(sparams, xx, [currV,currV,currV-0.2], 2, dotLocs, 0);
        
        diff = abs(currTc - tcBound);
    end
end

