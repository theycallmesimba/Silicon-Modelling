function [rhos, Hams] = simulateEffectiveShuttling(sparams, xx, vPulse, pTime, effHamParams, T2, detPulseFlag)
%SIMUALTEEFFECTIVESHUTTLING Summary of this function goes here
%   Detailed explanation goes here
    tic;
        
    % Are we given a voltage pulse or a detuning pulse as input
    if detPulseFlag
        epsL = vPulse(1,:);
        epsR = vPulse(2,:);
    else
        [epsL, epsR] = getDetuningVsVoltagePulse( sparams, xx, vPulse, 1 );
    end
    
%     figure;
%     hold on;
%     plot(epsL/sparams.ee,'linewidth',2);
%     plot(epsR/sparams.ee,'linewidth',2);
%     legend('epsL','epsR');
%     return
    
    % Testing Xuedong paper
%     pTime = 10E-9;
%     epsL = linspace(-0.2,0.2,length(vPulse(1,:)))*1E-3*sparams.ee;
%     epsR = linspace(0.2,-0.2,length(vPulse(1,:)))*1E-3*sparams.ee;
%     effHamParams{3} = 45E-6*sparams.ee;
    
    effHamParams{1} = epsL;
    effHamParams{2} = epsR;
    
    pulseTVec = linspace(0,pTime,length(vPulse(1,:)));
    
    detInterpolants = makeDetuningInterpolants(pulseTVec, epsL, epsR);

    time = 0:sparams.dt:pTime;
    fineEpsL = detInterpolants{1}({time});
    fineEpsR = detInterpolants{2}({time});

    h = waitbar(0,sprintf('Current Time Index: %d/%d',0,length(time)),...
        'Name','Simulating charge decoherence...',...
        'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');

    effHamParams{1} = fineEpsL(1);
    effHamParams{2} = fineEpsR(1);

    currRho = getEffectiveInitialRho(sparams,effHamParams);
    
    
    [rows,~] = size(currRho);
    if sparams.nStoreDataFrames > length(time)
        storeDataIndices = round(linspace(1,length(time),length(time)));
        rhos = zeros(rows,rows,length(time));
        Hams = zeros(rows,rows,length(time));
    else
        storeDataIndices = round(linspace(1,length(time),sparams.nStoreDataFrames));
        rhos = zeros(rows,rows,sparams.nStoreDataFrames);
        Hams = zeros(rows,rows,sparams.nStoreDataFrames);
    end
    
    ww = 1; % Store data indexing
    
    currHnp1 = constructEffectiveHamiltonian( sparams, effHamParams);
    
    rhos(:,:,ww) = currRho;
    Hams(:,:,ww) = currHnp1;
    for ii = 2:length(time)
        % Check for cancel button click
        if getappdata(h,'canceling')
            break;
        end

        % Update waitbar every N frames
        if mod(ii,5000) == 0
            waitbar(ii/length(time), h, sprintf('Current Time Index: %d/%d',ii,length(time)));
        end
    
        % Main part of simulation
        %*****************************************************************%
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
        %*****************************************************************%
        
        if any(storeDataIndices == ii)
            ww = ww + 1;
            rhos(:,:,ww) = currRho;
            
            Hams(:,:,ww) = currHn;
        end
    end
    
    delete(h);
    toc;
end

