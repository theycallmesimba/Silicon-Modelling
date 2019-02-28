function psi0 = getCoherentInitialState( sparams, xx, vPulse )
%GETINITIALSTATE Summary of this function goes here
%   Detailed explanation goes here

    % Solve the 1D SE for the initial potential well to get what our ground
    % state should look like
    vvInitial = sparams.P2DEGInterpolant([num2cell(vPulse(:,1)'),xx]);
    vvInitial = squeezeFast(length(sparams.gatesUsedInPulse),vvInitial);

    [psi0, ~] = solve1DSingleElectronSE(sparams,1,xx,vvInitial);
    
    if sparams.verbose
        fig = plotPotentialAndGroundWF(sparams,vPulse(:,1)',xx);
        title('Initial state','Fontsize',16);
        pause(1);
        delete(fig);
    end
end

