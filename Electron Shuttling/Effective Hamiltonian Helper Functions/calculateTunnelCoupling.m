function [tc, epsL, epsR] = calculateTunnelCoupling( sparams, xx, vv )
%CALCULATETUNNELCOUPLING Summary of this function goes here
%   Detailed explanation goes here
    epsL = 0;
    epsR = 0;

    % Find the ground state wave function for our potential
    [wf, ens] = solve1DSingleElectronSE(sparams, 2, xx, vv);
    dEstar = ens(2,2) - ens(1,1);
    
    % Now find out if the tunnel coupling is close to being turned on
    groundWF = wf(:,1).^2/getInnerProduct(xx,wf(:,1),wf(:,1));
    [wfPKS,~] = findpeaks(groundWF);
    wfPKS = wfPKS(wfPKS >= sparams.tcThreshold);
    if length(wfPKS) < 2
        [epsL, epsR] = getDetuning(sparams, xx, vv);
        tc = 0;
        % Now we need to find the 
        return;
    elseif length(wfPKS) > 2
        fprintf(1,'Found more than 2 peaks when trying to calculate tunnel coupling.\n');
        [epsL, epsR] = getDetuning(sparams, xx, vv);
        tc = 0;
        return;
    end
        
    [epsL, epsR] = getDetuning(sparams, xx, vv);
    
    dE = epsL - epsR;
    
    tc = sqrt(dEstar^2 - dE^2)/2;
    
    if ~isreal(tc)
        tc = 0;
    end
end

