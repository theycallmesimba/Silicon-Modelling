function [ sparams ] = getInitialState( sparams, xx )
%GETINITIALSTATE Summary of this function goes here
%   Detailed explanation goes here

    % Solve the 1D SE for the initial potential well to get what our ground
    % state should look like
    vvInitial = sparams.P2DEGInterpolant([num2cell(sparams.voltagePulse(:,1)'),xx]);
    vvInitial = squeezeFast(sparams.numOfGates,vvInitial);

    [sparams.rho0, ~] = solve1DSingleElectronSE(sparams,1,xx,vvInitial);
    
    if sparams.verbose
        fig = plotPotentialAndGroundWF(sparams,sparams.voltagePulse(:,1)',xx);
        title('Initial state','Fontsize',16);
        pause(5);
        delete(fig);
    end
end

