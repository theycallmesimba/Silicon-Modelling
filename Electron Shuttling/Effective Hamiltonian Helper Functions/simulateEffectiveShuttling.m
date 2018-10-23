function [fidelity, purity, orbExpValue, valExpValue, spinExpValue, finalRho] =...
    simulateEffectiveShuttling(sparams, xx, vPulse, pTime, effHamParams, T2, detPulseFlag)
%SIMUALTEEFFECTIVESHUTTLING Summary of this function goes here
%   Detailed explanation goes here
    tic;
    
    % Testing Xuedong paper
%     pTime = 10E-9;

    pulseTVec = linspace(0,pTime,length(vPulse(1,:)));

    % Are we given a voltage pulse or a detuning pulse as input
    if detPulseFlag
        epsL = vPulse(1,:);
        epsR = vPulse(2,:);
    else
        [epsL, epsR] = getDetuningVsVoltagePulse( sparams, xx, vPulse, 1 );
    end
    
    % Testing Xuedong paper
%     epsL = linspace(-0.2,0.2,length(vPulse(1,:)))*1E-3*sparams.ee;
%     epsR = linspace(0.2,-0.2,length(vPulse(1,:)))*1E-3*sparams.ee;
%     effHamParams{3} = 45E-6*sparams.ee;
    
    effHamParams{1} = epsL;
    effHamParams{2} = epsR;
    
    detInterpolants = makeDetuningInterpolants(pulseTVec, epsL, epsR);

    time = 0:sparams.dt:pTime;
    fineEpsL = detInterpolants{1}({time});
    fineEpsR = detInterpolants{2}({time});

    h = waitbar(0,sprintf('Current Time Index: %d/%d',0,length(time)),...
        'Name','Simulating charge decoherence...',...
        'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');

    effHamParams{1} = fineEpsL(1);
    effHamParams{2} = fineEpsR(1);
    
    H0 = constructEffectiveHamiltonian( sparams, effHamParams);
    [kets,ens] = eig(H0);
    [~, ind] = sort(diag(ens));
    kets = kets(:,ind);
    
    currRho0 = kets(:,1)*kets(:,1)';
    currRho = currRho0;
    
    % Testing Xuedong paper
%     kets(:,1) = kets(:,1)*exp(-1i*angle(kets(1,1)));
%     kets(:,2) = kets(:,2)*exp(-1i*angle(kets(1,2)));
%     currPsi = (kets(:,1) + kets(:,2))/sqrt(2);
%     currRho = currPsi*currPsi';
%     return
    
    fidelity = zeros(1,sparams.nFidelityFrames);
    purity = zeros(1,sparams.nPurityFrames);
    orbExpValue = zeros(3,sparams.nExpectationFrames);
    valExpValue = zeros(3,sparams.nExpectationFrames);
    spinExpValue = zeros(3,sparams.nExpectationFrames);
    
    fidIndices = round(linspace(1,length(time),sparams.nFidelityFrames));
    purIndices = round(linspace(1,length(time),sparams.nPurityFrames));
    expIndices = round(linspace(1,length(time),sparams.nExpectationFrames));
    
    nn = 1; % Fidelity indexing
    mm = 1; % Purity indexing
    ll = 1; % Orbital indexing
    fidelity(nn) = abs(trace(currRho0*currRho))^2; % Since we start in the ground state it should be 1
    purity(mm) = trace(currRho*currRho);
    [orbExp, valExp, spinExp] = getExpectationValues(sparams,currRho);
    orbExpValue(:,ll) = orbExp;
    valExpValue(:,ll) = valExp;
    spinExpValue(:,ll) = spinExp;
    
%     if sparams.includeValley && sparams.includeSpin
%         orbExpValue(1,ll) = trace(currRho*kron(sparams.sigmaz,eye(4)));
%         orbExpValue(2,ll) = trace(currRho*kron(sparams.sigmax,eye(4)));
%         orbExpValue(3,ll) = trace(currRho*kron(sparams.sigmay,eye(4)));
%     elseif sparams.includeValley || sparams.includeSpin
%         orbExpValue(1,ll) = trace(currRho*kron(sparams.sigmaz,eye(2)));
%         orbExpValue(2,ll) = trace(currRho*kron(sparams.sigmax,eye(2)));
%         orbExpValue(3,ll) = trace(currRho*kron(sparams.sigmay,eye(2)));
%     else
%         orbExpValue(1,ll) = trace(currRho*sparams.sigmaz);
%         orbExpValue(2,ll) = trace(currRho*sparams.sigmax);
%         orbExpValue(3,ll) = trace(currRho*sparams.sigmay);
%     end
    
    currHnp1 = H0;
    for ii = 2:length(time)
        % Check for cancel button click
        if getappdata(h,'canceling')
            break;
        end

        % Update waitbar every N frames
        if mod(ii,5000) == 0
            waitbar(ii/length(time), h, sprintf('Current Time Index: %d/%d',ii,length(time)));
        end

%         currH = updateEffectiveHamiltonianDetuning(currH,fineEpsL(ii-1)-fineEpsL(ii),fineEpsR(ii-1)-fineEpsR(ii));
        
%         currRho = currRho + sparams.dt*getMasterEquationOperator(sparams,currRho,currH,T2);
    
        % Attempt at Runge-Kutta method
        currHn = currHnp1;
        depsL = fineEpsL(ii) - fineEpsL(ii-1);
        depsR = fineEpsR(ii) - fineEpsR(ii-1);
        currHmid = updateEffectiveHamiltonianDetuning(currHn,depsL/2,depsR/2);
        currHnp1 = updateEffectiveHamiltonianDetuning(currHn,depsL,depsR);
        
        k1 = sparams.dt*getMasterEquationOperator(sparams,currRho,currHn,T2);
        k2 = sparams.dt*getMasterEquationOperator(sparams,currRho + k1/2,currHmid,T2);
        k3 = sparams.dt*getMasterEquationOperator(sparams,currRho + k2/2,currHmid,T2);
        k4 = sparams.dt*getMasterEquationOperator(sparams,currRho + k3,currHnp1,T2);
        
        currRho = currRho + 1/6*(k1 + 2*k2 + 2*k3 + k4);

        if any(fidIndices == ii)
            [kets,ens] = eig(currHnp1);
            [~, ind] = sort(diag(ens));
            kets = kets(:,ind);
            currRho0 = kets(:,1)*kets(:,1)';
            
            nn = nn + 1;
            fidelity(nn) = abs(trace(currRho0*currRho))^2;
        end
        
        if any(purIndices == ii)
            mm = mm + 1;
            purity(mm) = trace(currRho*currRho);
        end
        
        if any(expIndices == ii)
            ll = ll + 1;
            [orbExp, valExp, spinExp] = getExpectationValues(sparams,currRho);
            orbExpValue(:,ll) = orbExp;
            valExpValue(:,ll) = valExp;
            spinExpValue(:,ll) = spinExp;
        end
    end
    
    finalRho = currRho;

    delete(h);
    toc;
end

