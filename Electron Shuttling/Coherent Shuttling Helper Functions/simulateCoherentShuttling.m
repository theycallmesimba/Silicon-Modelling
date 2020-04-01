function sparams = simulateCoherentShuttling(sparams, sweepVec, xx, zz, pp, currSimIndex, simMethod)
%SIMULATECOHERENTSHUTTLING Summary of this function goes here
%   Detailed explanation goes here

    % Make folder to save images for individual simulations
    if strcmp(sparams.sweptParameter,'time')        
        currSimFolder = num2str(sweepVec(currSimIndex));
    elseif strcmp(sparams.sweptParameter,'adiabicity')     
        currSimFolder = [num2str(sweepVec(currSimIndex,1)) '-' num2str(sweepVec(currSimIndex,2))];
    end
    mkdir([sparams.saveDir sparams.saveFolder currSimFolder]);

    if strcmp(simMethod,'split-operator')
        [~, ~, K, K2] = initializeShuttlingSimulation(sparams, pp);
    end
    
    % Now we need to make the individual gate interpolants for the pulse
    % First, we want to associate each potential simulation we have with a time
    % value (i.e. when in the simulation should that potential appear)
    tPots = linspace(0,sparams.totalTime(currSimIndex),length(sparams.voltagePulse(currSimIndex,1,:)));
    tTime = 0:sparams.dt:sparams.totalTime(currSimIndex);

    % Build the pulse interpolants for this simulation
    sparams.vPulseGInterpolants = makePulseInterpolants(sparams, tPots,...
        squeeze(sparams.voltagePulse(currSimIndex,:,:)));

    % Get time indices to save figures
    sparams.saveFigureIndices(currSimIndex,:) = round(linspace(1,length(tTime),sparams.nFigureFrames));
    % Get time indices and corresponding time values to calculate and save starkShift
    sparams.starkShiftIndices(currSimIndex,:) = round(linspace(1,length(tTime),sparams.nStarkShiftFrames));
    sparams.tStarkShift(currSimIndex,:) = tTime(sparams.starkShiftIndices(currSimIndex,:));
    sparams.fidelityIndices(currSimIndex,:) = round(linspace(1,length(tTime),sparams.nStoreDataFrames));
    
    % Make waitbar
    if strcmp(sparams.sweptParameter,'time')        
        % Make the waitbar to show run time
        h = waitbar(0,sprintf('Current Time Index: %d/%d',0,length(tTime)),...
            'Name',sprintf('Performing time shuttling simulation for %E...',sweepVec(currSimIndex)),...
            'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
    elseif strcmp(sparams.sweptParameter,'adiabicity')        
        % Make the waitbar to show run time
        h = waitbar(0,sprintf('Current Time Index: %d/%d',0,length(tTime)),...
            'Name',sprintf('Performing adiabatic shuttling simulation for %E...',sweepVec(currSimIndex)),...
            'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
    end
    movegui(h,'northwest');
    
    % Set current Psi to be the initial state
    currPsi = getCoherentInitialState(sparams,xx,squeeze(sparams.voltagePulse(currSimIndex,:,:)));
            
    % Initialize the shuttling simulation figure
    shtlEvolutionFig = initializeShuttlingFigure(sparams, squeeze(sparams.voltagePulse(currSimIndex,:,:)),...
        currPsi, currPsi, xx, sparams.totalTime(currSimIndex));    
        
    nn = 1; % Used to index fidelity array
    ll = 0; % Used to know where in time domain to interpolate our potentials
    kk = 0; % Used to index which interpolated potential we are on
    yy = 0; % Used for stark shift indexing

    switch simMethod
        case 'split-operator'
            % Convert from position to momentum space
            currPsip = fftshift(fft(currPsi));
            % Apply the KE operator for dt/2 
            currPsip = K.*currPsip; 
    end
    
    for ii = 1:length(tTime)
        kk = kk + 1;
        
        % Check for cancel button click
        if getappdata(h,'canceling')
            flag = 1;
            break;
        end
        
        % Update waitbar every N frames
        if mod(ii,sparams.updateWaitbar) == 0
            waitbar(ii/length(tTime), h, sprintf('Current Time Index: %d/%d',ii,length(tTime)));
        end

        % Get an updated set of voltage pulse values 
        if mod(ii,sparams.updateInterpPot) == 0 || ii == 1
            kk = 1; % Reset counter
            
            startInterpInd = ll*sparams.updateInterpPot;
            if ii == 1
                startInterpInd = 1;
            end
            
            % Increment counter for what batch of time indices in the pulse
            % to interpolate
            ll = ll + 1;
            
            endInterpInd = ll*sparams.updateInterpPot - 1;
            if endInterpInd > length(tTime)
                endInterpInd = length(tTime);
            end
            
            gPulse = getInterpolatedPulseValues(sparams,...
                tTime(startInterpInd:endInterpInd),sparams.vPulseGInterpolants);
        end
        
        currPotential = sparams.P2DEGInterpolant(getInterpolantArgument(gPulse(:,kk),xx));
        currPotential = squeezeFast(length(sparams.gatesUsedInPulse),currPotential)';
        
        % Main part of simulation
        %*****************************************************************%
        switch simMethod
            case 'split-operator'
                V = exp(-1i*sparams.dt*currPotential'/sparams.hbar);

                % Convert from momentum to position space
                currPsix = ifft(fftshift(currPsip));
                % Apply the PE operator for dt
                currPsix = V.*currPsix;
                % Convert from position to momentum space
                currPsip = fftshift(fft(currPsix));
                if ii ~= length(tTime)
                    % Apply the KE operator for dt
                    currPsip = K2.*currPsip;
                else
                    % Apply the KE operator for dt/2
                    currPsip = K.*currPsip;
                    % Convert from momentum to position space
                    currPsi = ifft(fftshift(currPsip));
                end
            case 'runge-kutta'
        end
        %*****************************************************************%

        % Calculate Stark shift
        if any(sparams.starkShiftIndices(currSimIndex,:) == ii) && sparams.calculateStarkShift
            yy = yy + 1;

            curr2DPot = squeezeFast(length(sparams.gatesUsedInPulse),...
                sparams.P2DInterpolant(getInterpolantArgument(gPulse(:,kk),xx,zz)));
            switch simMethod
                case 'split-operator'
                    sparams = calculateStarkShift(sparams,curr2DPot,...
                        ifft(fftshift(currPsip)),xx,currSimIndex,yy);
                case 'runge-kutta'
                    sparams = calculateStarkShift(sparams,curr2DPot,...
                        currPsi,xx,currSimIndex,yy);
            end
            sparams = calculateStarkShift(sparams,curr2DPot,...
                ifft(fftshift(currPsip)),xx,currSimIndex,yy);
        end
        
        % Update figure 
        if mod(ii,sparams.updateFigure) == 0
            [currPsi0, ~] = solve1DSingleElectronSE(sparams,1,xx,currPotential);
            
            switch simMethod
                case 'split-operator'
                    updateShuttlingFigure(sparams,shtlEvolutionFig,ifft(fftshift(currPsip)),...
                        currPsi0,currPotential);
                case 'runge-kutta'
                    updateShuttlingFigure(sparams,shtlEvolutionFig,currPsi,...
                        currPsi0,currPotential);
            end
        end
        
        % Update figure and save to gif
        if any(sparams.saveFigureIndices(currSimIndex,:) == ii)
            [currPsi0, ~] = solve1DSingleElectronSE(sparams,1,xx,currPotential);

            switch simMethod
                case 'split-operator'
                    updateShuttlingFigure(sparams,shtlEvolutionFig,ifft(fftshift(currPsip)),...
                        currPsi0,currPotential);
                case 'runge-kutta'
                    updateShuttlingFigure(sparams,shtlEvolutionFig,currPsi,...
                        currPsi0,currPotential);
            end
            
            saveGIFofEvolution(shtlEvolutionFig, sweepVec(currSimIndex), tTime(ii),...
                [sparams.saveDir sparams.saveFolder currSimFolder]);
        end
        
        % Calculate fidelity WRT current ground state
        if any(sparams.fidelityIndices(currSimIndex,:) == ii)            
            % Need to get the ground state of the current potential
            [currPsi0, ~] = solve1DSingleElectronSE(sparams,1,xx,currPotential);
            
            switch simMethod
                case 'split-operator'
                    sparams.fidelity(currSimIndex,nn) = abs(getInnerProduct(xx,currPsi0,ifft(fftshift(currPsip)))).^2;
                case 'runge-kutta'
                    sparams.fidelity(currSimIndex,nn) = abs(getInnerProduct(xx,currPsi0,currPsi)).^2;
            end
            nn = nn + 1;
        end
    end
        
    % Close simulation figure
    close(shtlEvolutionFig);
    % Close waitbar
    delete(h);
end

