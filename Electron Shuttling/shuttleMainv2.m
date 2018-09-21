% Load all the parameters for the simulation
clear sparams xx vv pp;
shuttleParameterFile;

fprintf(1,'Loading potentials...\n');
[sparams,xx,zz] = loadPotentials(sparams);

sparams.nxGrid = length(xx);
sparams.nzGrid = length(zz);
sparams.dx = xx(2) - xx(1);
sparams.dz = zz(2) - zz(1);
sparams.dp = 2*pi*sparams.hbar/(sparams.dx*sparams.nxGrid);
pp = ((-sparams.nxGrid/2):1:(sparams.nxGrid/2 - 1))*sparams.dp;

% Find which index corresponds to where the 2DEG should be
[~,sparams.twoDEGindZ] = min(abs(zz - (-0.5*1E-9)));
for ii = 1:length(sparams.potentials)
    sparams.potentials(ii).pot2DEG = sparams.potentials(ii).pot2D(sparams.twoDEGindZ,:);
end

% Now we want to make the potential interpolant object (both 2D and 2DEG)
% sparams = makePotentialsInterpolants(sparams,xx,zz);
sparams = makePotentialsInterpolants(sparams,xx,zz);

% Check that the potentials were loaded correctly and that the interpolants
% were correctly assembled
debugHere = 0;
if debugHere
    checkPotentialLoad(sparams,xx,zz);
end
%%
% Get our desired votlage pulse.
% vec = [[0.01,0.01];[0.001,0.025];[0.001,0.05]];
vec = [0.05,0.05];
for ii = 1:length(vec(:,1))
    [sparams, vPulse, vPulseTime] = getVoltagePulseAdiabatic(sparams,xx,vec(ii,:),0.8);
    vPulseTime
end
sparams.voltagePulse = vPulse;

debugHere = 0;
if debugHere
    analyzePulseAdiabicity(sparams,xx,squeeze(vPulse(1,:,:)),...
        vPulseTime(1),1,1:4);
end
%%
profile on;

% Using ref https://arxiv.org/pdf/1306.3247.pdf we now find the time
% evolution operator U(t + dT,t) to evolve our initial wavefunction to the
% next wavefunction.  That newly found wavefunction will act as our initial
% state for the next time frame.  This method uses the split ooperator
% approach
[sparams, sweepVec, K, K2] = initializeShuttlingSimulation(sparams, pp);

for jj = 1:length(sweepVec)
    if strcmp(sparams.sweptParameter,'time')
        fprintf(1,'Running time sweep shuttling simulation for %.3E (%d/%d)...\n',sweepVec(jj),jj,length(sweepVec));
    elseif strcmp(sparams.sweptParameter,'adiabicity')
        fprintf(1,'Running adiabatic sweep shuttling simulation for [%.3E,%.3E] (%d/%d)...\n',sweepVec(jj,1),sweepVec(jj,2),jj,length(sweepVec));
    end
    
    % Make folder to save images for individual simulations
    if strcmp(sparams.sweptParameter,'time')        
        currSimFolder = num2str(sweepVec(jj));
    elseif strcmp(sparams.sweptParameter,'adiabicity')     
        currSimFolder = [num2str(sweepVec(jj,1)) '-' num2str(sweepVec(jj,2))];
    end
    mkdir([sparams.saveDir sparams.saveFolder currSimFolder]);
    
    tic;
    
    fprintf(1,'Building voltage pulse...\n');
    % First step is to get the voltage pulse depending on our simulation
    % type
    if strcmp(sparams.sweptParameter,'time') 
        [sparams, sparams.voltagePulse(jj,:,:)] = getVoltagePulse(sparams,xx);
    elseif strcmp(sparams.sweptParameter,'adiabicity')
        [sparams, temp, sparams.totalTime(jj)] = getVoltagePulseAdiabatic(sparams,xx,sweepVec(jj,:));
        sparams.voltagePulse(jj,:,:) = temp;
    end
    
    fprintf(1,'Getting initial wavefunction...\n');
    sparams = getInitialState(sparams,xx,squeeze(sparams.voltagePulse(jj,:,:)));
    
    % Now we need to make the individual gate interpolants for the pulse
    % First, we want to associate each potential simulation we have with a time
    % value (i.e. when in the simulation should that potential appear)
    tPots = linspace(0,sparams.totalTime(jj),length(sparams.voltagePulse(jj,1,:)));
    tTime = 0:sparams.dt:sparams.totalTime(jj);

    % Build the pulse interpolants for this simulation
    sparams.vPulseGInterpolants = makePulseInterpolants(sparams, tPots,...
        squeeze(sparams.voltagePulse(jj,:,:)));

    % Get time indices to save figures
    sparams.saveFigureIndices(jj,:) = round(linspace(1,length(tTime),sparams.nFigureFrames));
    % Get time indices and corresponding time values to calculate and save starkShift
    sparams.starkShiftIndices(jj,:) = round(linspace(1,length(tTime),sparams.nStarkShiftFrames));
    sparams.tStarkShift(jj,:) = tTime(sparams.starkShiftIndices(jj,:));
    sparams.fidelityIndices(jj,:) = round(linspace(1,length(tTime),sparams.nFidelityFrames));
    
    % Make waitbar
    if strcmp(sparams.sweptParameter,'time')        
        % Make the waitbar to show run time
        h = waitbar(0,sprintf('Current Time Index: %d/%d',0,length(tTime)),...
            'Name',sprintf('Performing time shuttling simulation for %E...',sweepVec(jj)),...
            'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
    elseif strcmp(sparams.sweptParameter,'adiabicity')        
        % Make the waitbar to show run time
        h = waitbar(0,sprintf('Current Time Index: %d/%d',0,length(tTime)),...
            'Name',sprintf('Performing adiabatic shuttling simulation for %E...',sweepVec(jj)),...
            'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
    end
    movegui(h,'northwest');

    % Set current Psi to be the initial state
    currPsi = sparams.rho0';
    
    % Initialize the shuttling simulation figure
    shtlEvolutionFig = initializeShuttlingFigure(sparams, squeeze(sparams.voltagePulse(jj,:,:)),...
        currPsi, currPsi, xx, sparams.totalTime(jj));
    
    nn = 1; % Used to index fidelity array
    ll = 0; % Used to know where in time domain to interpolate our potentials
    kk = 0; % Used to index which interpolated potential we are on
    yy = 0; % Used for stark shift indexing

    % Convert from position to momentum space
    currPsip = fftshift(fft(currPsi));
    % Apply the KE operator for dt/2 
    currPsip = K.*currPsip;
    
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
        currPotential = squeezeFast(sparams.numOfGates,currPotential)';
        V = exp(-1i*sparams.dt*currPotential/sparams.hbar);
        
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
        
        % Calculate Stark shift
        if any(sparams.starkShiftIndices(jj,:) == ii) && sparams.calculateStarkShift
            yy = yy + 1;

            curr2DPot = squeezeFast(sparams.numOfGates,...
                sparams.P2DInterpolant(getInterpolantArgument(gPulse(:,kk),xx,zz)));
            sparams = calculateStarkShift(sparams,curr2DPot,...
                ifft(fftshift(currPsip)),xx,jj,yy);
        end
        
        % Update figure 
        if mod(ii,sparams.updateFigure) == 0
            [currRho0, ~] = solve1DSingleElectronSE(sparams,1,xx,currPotential);
            
            updateShuttlingFigure(sparams,shtlEvolutionFig,ifft(fftshift(currPsip)),...
                currRho0,currPotential);
        end
        
        % Update figure and save to gif
        if any(sparams.saveFigureIndices(jj,:) == ii)
            [currRho0, ~] = solve1DSingleElectronSE(sparams,1,xx,currPotential);
            
            updateShuttlingFigure(sparams,shtlEvolutionFig,ifft(fftshift(currPsip)),...
                currRho0,currPotential);
            
            saveGIFofEvolution(shtlEvolutionFig, sweepVec(jj), tTime(ii),...
                [sparams.saveDir sparams.saveFolder currSimFolder]);
        end
        
        % Calculate fidelity WRT current ground state
        if any(sparams.fidelityIndices(jj,:) == ii)            
            % Need to get the ground state of the current potential
            [currRho0, ~] = solve1DSingleElectronSE(sparams,1,xx,currPotential);
            sparams.fidelity(jj,nn) = abs(getInnerProduct(xx,currRho0.',ifft(fftshift(currPsip))))^2;
            nn = nn + 1;
        end
        
    end
            
    % Close simulation figure
    close(shtlEvolutionFig);
    % Close waitbar
    delete(h);
    
    % Save the simulation results so far (in case of a shut down mid
    % simulation)
    save([sparams.saveDir sparams.saveFolder 'simulationData'],'sparams','xx','zz','pp');
    
    toc;
end

profile off
profile viewer
%%
% Post simulation Analysis
analyzeFidelity(sparams)
analyzeStarkShift(sparams)
%% Post simulation Analysis

figure;
plot(linspace(0,4,sparams.nStarkShiftFrames),(sparams.vShiftGround(23,:) - sparams.vShiftGround(23,1))*1E-6,'Linewidth',2.5);
line([0,4],[mean(sparams.vShiftGround(23,:) - sparams.vShiftGround(23,1)),mean(sparams.vShiftGround(23,:) - sparams.vShiftGround(23,1))]*1E-6,...
    'Linewidth',2.5,'Linestyle','--','Color','red');
set(gca,'Fontsize',15);
xlabel('Time [ns]','Interpreter','Latex','Fontsize',24);
ylabel('$\nu - \nu_0$ [MHz]','Interpreter','Latex','Fontsize',24);

%%
% Need to get the Singlet decay
% Get singlet fidelity decay

rho0 = [0;1;-1;0]/sqrt(2);
rhoIdeal = rho0;
rhoStark = rho0;
rhoIdealMean = rho0;
rhoStarkMean = rho0;
tIndex = 13;

times = linspace(0,4E-9,sparams.nStarkShiftFrames + 1);
dt = times(2) - times(1);
dz = [1,0;0,-1];
dz1 = kron(dz,eye(2));
dz2 = kron(eye(2),dz);
v1 = sparams.vShiftGround(tIndex,1);
v0mean = mean(sparams.vShiftGround(tIndex,:));
fidelity(1) = 1;
fidelityMean(1) = 1;
for ii = 1:length(times)-1
    Uideal = expm(-1i*dt*pi*(dz1*v1 + dz2*v1));
    Ustark = expm(-1i*dt*pi*(dz1*v1 + dz2*sparams.vShiftGround(tIndex,ii)));
    UidealMean = expm(-1i*dt*pi*(dz1*v0mean + dz2*v0mean));
    UstarkMean = expm(-1i*dt*pi*(dz1*v0mean + dz2*sparams.vShiftGround(tIndex,ii)));

    rhoIdeal = Uideal*rhoIdeal;
    rhoStark = Ustark*rhoStark;
    
    rhoIdealMean = UidealMean*rhoIdealMean;
    rhoStarkMean = UstarkMean*rhoStarkMean;
    
    fidelity(ii+1) = abs(rhoIdeal'*rhoStark)^2;
    fidelityMean(ii+1) = abs(rhoIdealMean'*rhoStarkMean)^2;
end

figure;
hold on;
plot(times*1E9,1-fidelity,'Linewidth',2.5);
plot(times*1E9,1-fidelityMean,'Linewidth',2.5,'Color','red');
% set(gca, 'YTick', [10^-15, 10^-13, 10^-11, 10^-9, 10^-7]);
% yticks(logspace(-15,-9,7));
set(gca, 'YScale', 'log');
grid minor;
grid;
set(gca,'Fontsize',14);
xlabel('Time [ns]','Interpreter','Latex','Fontsize',22);
ylabel('$1 - |\langle S|\psi(t)\rangle|^2$','Interpreter','Latex','Fontsize',22);
% set(gca,'Fontsize',14);
% ax = gca;
% set(gca,'YTickLabel',sprintf('%.2e\n',fidelity))
%%
figure;
hold on;
set(gca,'Fontsize',14);
for ii = 1:sparams.numOfGates
    plot(sparams.voltagePulse(ii,:),'Linewidth',2);
    xlabel('Percentage of total shuttling time [\%]','Interpreter','Latex','Fontsize',22);
    ylabel('Voltage [V]','Interpreter','Latex','Fontsize',22);
    xlim([1,length(sparams.voltagePulse(1,:))]);
    ylim([min(min(sparams.voltagePulse)),max(max(sparams.voltagePulse))*1.01]);
end
legend(sparams.gateLabels);

timeIndex = 3;
gateVoltages = [sparams.voltagePulse(1,timeIndex),...
    sparams.voltagePulse(2,timeIndex), sparams.voltagePulse(3,timeIndex)];
tempPot = squeeze(sparams.P2DEGInterpolant([num2cell(gateVoltages),...
        mat2cell(xx,1,length(xx))]));
[tempWF, ~] = solve1DSingleElectronSE(sparams,1,xx,tempPot); 
tempWFNSqr = abs(tempWF).^2/norm(abs(tempWF).^2);

% figure;
% plot(xx,tempWFNSqr);
% pks = findpeaks(tempWFNSqr);
% pks


%%
thresh = linspace(0.00005,0.0015,15);
plot(thresh,temp,'Linewidth',1.6);
set(gca,'Fontsize',13);
xlabel('Peak Height Threshold','Interpreter','Latex','Fontsize',22);
ylabel('Voltage [V]','Interpreter','Latex','Fontsize',22);

%%
sweepVec = getSweepVector(sparams);
figure;
plot(sweepVec,sparams.tStarkShift(:,500));
    
    
    