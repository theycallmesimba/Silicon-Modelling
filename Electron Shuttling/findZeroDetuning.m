function vGTarg = findZeroDetuning( sparams, xx, vVec, gateIndSweep )
%FINDZERODETUNING Finds what gate voltage value will give you zero detuning
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
    options = optimset('TolX',1E-6,'Display','iter');
    [vGTarg, ~] = fminbnd(@(x) findMinDeltaE(x),minGVal,maxGVal,options);
    
    function deltaE = findMinDeltaE(gateValue)
        % Update the current voltage vector and get the corresponding
        % potential
        vVecCell{gateIndSweep} = gateValue;
        currPot = squeeze(sparams.P2DEGInterpolant(vVecCell));
        
        % Find the wells of the potential calculate their detuning
        peaks = sort(findpeaks(-currPot),'descend');
        deltaE = (peaks(1) - peaks(2))/sparams.ee;
    end
end

