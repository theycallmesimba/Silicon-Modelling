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
    
    % Now, the SE equation solver uses eigs which can output a random +/-1
    % global phase on the wavefunction output depending on what random seed
    % vector it uses to solve the eigenvalue problem.  We need to make sure
    % our phases are consistent across all dts we calculate otherwise we'll
    % get big jumps in the derivative.  We can check the phase by taking
    % the inner product which (if dt is small) should be ~1.  So if it's <0
    % then we know the phase is off.
    [~,nBasis] = size(currHam);
    for ii = 1:nBasis
        if kets(:,ii)'*ketsp2h(:,ii) < 0
            ketsp2h(:,ii) = -ketsp2h(:,ii);
        end
        if kets(:,ii)'*ketsph(:,ii) < 0
            ketsph(:,ii) = -ketsph(:,ii);
        end
        if kets(:,ii)'*ketsmh(:,ii) < 0
            ketsmh(:,ii) = -ketsmh(:,ii);
        end
        if kets(:,ii)'*ketsm2h(:,ii) < 0
            ketsm2h(:,ii) = -ketsm2h(:,ii);
        end
    end
    
    
    % Now evaluate the adiabatic parameter by summing over all excited
    % states
    adiabParamNNTotal = 0;
    vv = 0;
    adiabParamInd = struct([]);
    % <jj|\dot{ii}> adiabatic parameter
    for ii = 1:nBasis
        % Get |\dot{ii}>
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
            
            if ii == nn
                adiabParamNNTotal = adiabParamNNTotal + adiabParamInd(vv).adiabaticParameter;
            end
        end
    end
end