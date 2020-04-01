function [adiabParamNNTotal, adiabParamInd] = calculateAdiabaticParameterEffective(...
    sparams, detInterpolants, cTime, nn, effHamiltonianParams )
    
    h = sparams.hdt;

    % Get current pulse point plus the shifted ones as well (+h,+2h,-h,-2h)
    [epsLPoints, epsRPoints] = getInterpolatedDetuningValues(...
        [cTime-(2*h),cTime-h,cTime,cTime+h,cTime+(2*h)],detInterpolants);
    
    % Solve the Hamiltonian at the current time point to get vectors and 
    effHamiltonianParams{1} = epsLPoints(3);
    effHamiltonianParams{2} = epsRPoints(3);
    currHam = constructEffectiveHamiltonian( sparams, effHamiltonianParams);
    [kets,ens] = eig(currHam);
    [~, ind] = sort(diag(ens));
    ens = ens(ind,ind);
    kets = kets(:,ind);
    
    % Estimate the derivative of the nth WF using a five point stencil
    effHamiltonianParams{1} = epsLPoints(5);
    effHamiltonianParams{2} = epsRPoints(5);
    currHamp2h = constructEffectiveHamiltonian( sparams, effHamiltonianParams);
    [ketsp2h,ensp2h] = eig(currHamp2h);
    [~, ind] = sort(diag(ensp2h));
    ketsp2h = ketsp2h(:,ind);
    
    effHamiltonianParams{1} = epsLPoints(4);
    effHamiltonianParams{2} = epsRPoints(4);
    currHamph = constructEffectiveHamiltonian( sparams, effHamiltonianParams);
    [ketsph,ensph] = eig(currHamph);
    [~, ind] = sort(diag(ensph));
    ketsph = ketsph(:,ind);
    
    effHamiltonianParams{1} = epsLPoints(2);
    effHamiltonianParams{2} = epsRPoints(2);
    currHammh = constructEffectiveHamiltonian( sparams, effHamiltonianParams);
    [ketsmh,ensmh] = eig(currHammh);
    [~, ind] = sort(diag(ensmh));
    ketsmh = ketsmh(:,ind);
    
    effHamiltonianParams{1} = epsLPoints(1);
    effHamiltonianParams{2} = epsRPoints(1);
    currHamm2h = constructEffectiveHamiltonian( sparams, effHamiltonianParams);
    [ketsm2h,ensm2h] = eig(currHamm2h);
    [~, ind] = sort(diag(ensm2h));
    ketsm2h = ketsm2h(:,ind);
    
    % Now, the SE equation solver uses eigs which can output a random
    % global phase on the eigenvector outputs.  This will cause issues when
    % we try to find the derivative as they need the same global phase in
    % order to accurately add them together.  So we normalize them all to
    % the same phase by making the first element of every vector real and
    % positive
    ketsp2h = setKetSameGlobalPhase(ketsp2h);
    ketsph = setKetSameGlobalPhase(ketsph);
    ketsmh = setKetSameGlobalPhase(ketsmh);
    ketsm2h = setKetSameGlobalPhase(ketsm2h);
    
    [~,nBasis] = size(currHam);
    % Now evaluate the adiabatic parameter by summing over all excited
    % states
    adiabParamNNTotal = 0;
    vv = 0;
    adiabParamInd = struct([]);
    % <jj|\dot{ii}> adiabatic parameter
    for ii = 1:nBasis
        % Get |\dot{ii}> through a 5 point stencil estimation of the
        % derivative
        iiKetDeriv = (-ketsp2h(:,ii) + 8*ketsph(:,ii) - 8*ketsmh(:,ii) + ketsm2h(:,ii))/(12*h);
        for jj = 1:nBasis
            % Ignore <ii|\dot{ii}> terms
            if ii == jj
                continue
            end
            vv = vv + 1;
            adiabParamInd(vv).mm = jj;
            adiabParamInd(vv).nn = ii;
            adiabParamInd(vv).adiabaticParameter = sparams.hbar*abs(kets(:,jj)'*iiKetDeriv/(ens(ii,ii) - ens(jj,jj)));
            
            if any(ii == nn)
                adiabParamNNTotal = adiabParamNNTotal + adiabParamInd(vv).adiabaticParameter;
            end
        end
    end
    
    % Now, we want to make sure that a certain adiabatic parameter
    % corresponds to a certain fidelity..  So we will average the adiabatic
    % parameter to how many nn indices we swept over.
    adiabParamNNTotal = adiabParamNNTotal/length(nn);
end