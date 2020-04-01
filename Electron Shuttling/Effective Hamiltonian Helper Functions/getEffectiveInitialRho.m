function rho0 = getEffectiveInitialRho(sparams,effHamParams)
%GETEFFECTIVEINITIALSTATE Summary of this function goes here
%   Detailed explanation goes here

    if sparams.includeSecondSpin
        % First, just get ground state of valley and orbit components
        sparams.includeSpin = 0;
        sparams.includeSecondSpin = 0;
        H0 = constructEffectiveHamiltonian( sparams, effHamParams);
        [kets,ens] = eig(H0);
        [~, ind] = sort(diag(ens));
        kets = kets(:,ind);
        kets = setKetSameGlobalPhase(kets);
        
        psiOrbVal0 = kets(:,1);
        rhoOrbVal0 = psiOrbVal0*psiOrbVal0';
        
        % DELETE THIS SECTION LATER
        sparams.includeSpin = 1;
        sparams.includeOrbital = 0;
        sparams.includeValley = 0;
        sparams.includeSecondSpin = 1;
        
        
        psiSinglet = 1/sqrt(2)*[0,1,-1,0]';
        rhoSinglet = psiSinglet*psiSinglet';
        
        rho0 = kron(rhoOrbVal0,rhoSinglet);
        return;
    end
    
    sparams.includeSpin = 1;
    sparams.includeSecondSpin = 0;
    H0 = constructEffectiveHamiltonian( sparams, effHamParams);
    [kets,ens] = eig(H0);
    [~, ind] = sort(diag(ens));
    kets = kets(:,ind);
    kets = setKetSameGlobalPhase(kets);
    
    [rows,~] = size(kets);
    psi0 = zeros(rows,1);
    for ii = 1:length(sparams.stateIndices)
        psi0 = psi0 + 1/sqrt(length(sparams.stateIndices))*kets(:,ii);
    end
    
    rho0 = psi0*psi0';
end

