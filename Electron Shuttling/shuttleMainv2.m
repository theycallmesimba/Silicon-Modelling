% Load all the parameters for the simulation
clear sparams xx vv;
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
% Plot a random potential and the ground state
debugHere = 1;
if debugHere
    % Plot a random potential and the ground state
    plotPotentialAndGroundWF(sparams,[0.8,0.7998,0.6,0.6,0.8],xx);
end

%%
% Get our desired votlage pulse.
sparams = getVoltagePulse(sparams,xx);

debugHere = 1;
if debugHere
    fig = checkVoltagePulse(sparams);
    pause(5);
    delete(fig);
end
%%
% analyzeEEnergiesVersusPulse(sparams, xx);

profile on;

fprintf(1,'Getting initial wavefunction...\n');
sparams = getInitialState(sparams,xx);


% Using ref https://arxiv.org/pdf/1306.3247.pdf we now find the time
% evolution operator U(t + dT,t) to evolve our initial wavefunction to the
% next wavefunction.  That newly found wavefunction will act as our initial
% state for the next time frame.  This method uses the split ooperator
% approach
    
% Make the KE operator since it is the same every time it is applied
K = exp(-1i*sparams.dt/2*(pp.^2)/(2*sparams.me*sparams.hbar));
K2 = K.*K;

% Make the fidelity array
% maxTime = max(sparams.totalTime);
% maxLength = length(0:sparams.dt:maxTime);
sparams.fidelity = zeros(length(sparams.totalTime),sparams.nFidelityFrames);
sparams.starkShift = zeros(length(sparams.totalTime),sparams.nStarkShiftFrames);

sparams.avgEzGround = zeros(length(sparams.totalTime),sparams.nStarkShiftFrames);
sparams.avgEz = zeros(length(sparams.totalTime),sparams.nStarkShiftFrames);
sparams.vShiftGround = zeros(length(sparams.totalTime),sparams.nStarkShiftFrames);
sparams.vShift = zeros(length(sparams.totalTime),sparams.nStarkShiftFrames);

% Let's create the folder to save data
time = clock;
sparams.saveFolder = sprintf('%d-%02d-%02d-%02d-%02d-%02d',time(1),time(2),time(3),time(4),time(5),round(time(6)));
mkdir([sparams.saveDir sparams.saveFolder]);
sparams.saveFolder = [sparams.saveFolder '/'];

for jj = 1:length(sparams.totalTime)
    mkdir([sparams.saveDir sparams.saveFolder num2str(sparams.totalTime(jj))]);
    
    tic;
    
    % Now we need to make the individual gate interpolants for the pulse
    % First, we want to associate each potential simulation we have with a time
    % value (i.e. when in the simulation should that potential appear)
    tPots = linspace(0,sparams.totalTime(jj),length(sparams.voltagePulse(1,:)));
    sparams.vPulseGInterpolants = {};
    for vv = 1:sparams.numOfGates
        sparams.vPulseGInterpolants{vv} = griddedInterpolant({tPots},sparams.voltagePulse(vv,:));
    end
    
    % Get number of time steps
    tTime = 0:sparams.dt:sparams.totalTime(jj);

    % Get time indices to save figures
    sparams.saveFigureIndices(jj,:) = round(linspace(1,length(tTime),sparams.nFigureFrames));
    % Get time indices and corresponding time values to calculate and save starkShift
    sparams.starkShiftIndices(jj,:) = round(linspace(1,length(tTime),sparams.nStarkShiftFrames));
    sparams.tStarkShift(jj,:) = tTime(sparams.starkShiftIndices(jj,:));
    sparams.fidelityIndices(jj,:) = round(linspace(1,length(tTime),sparams.nFidelityFrames));
    
    fprintf(1,'Running shuttling simulation for %E (%d/%d)...\n',sparams.totalTime(jj),jj,length(sparams.totalTime));

    % Make the waitbar to show run time
    h = waitbar(0,sprintf('Current Time Index: %d/%d',0,length(tTime)),...
        'Name',sprintf('Performing shuttling simulation for %E...',sparams.totalTime(jj)),...
        'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
    movegui(h,'northwest');

    currPsi = sparams.rho0';
    
    shtlEvolutionFig = initializeShuttlingFigure(sparams, currPsi, currPsi, xx, jj);
    
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
                        
            saveGIFofEvolution(sparams, shtlEvolutionFig, tTime(ii), sparams.totalTime(jj));
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
    save([sparams.saveDir sparams.saveFolder 'simulationData'],'sparams','xx','zz');
    
    toc;
end

profile off
profile viewer

%% Post simulation Analysis
% fids = sparams.fidelity;
fids = sparams.fidelity;
fids(fids==0) = NaN;
[rows,cols] = size(fids);
highTime = max(sparams.totalTime);
fidTimeIndices = sparams.updateFidelity*sparams.dt:sparams.updateFidelity*sparams.dt:highTime;
[TIndex,TTime] = meshgrid(fidTimeIndices,[0,sparams.totalTime,2*max(sparams.totalTime)]);
fidelTemp = zeros(rows+2,cols);
fidelTemp(2:(rows+1),:) = fids;

figure;
s = surf(TIndex,TTime,fidelTemp);
set(s,'edgecolor','none');
set(gca,'XScale','log');
set(gca,'YScale','log');
xlabel('Time step $t_j$ [s]','interpreter','latex','fontsize',15);
ylabel('Total Shuttling Simulated Time [s]','interpreter','latex','fontsize',15);
xlim([min(min(TIndex)),max(max(TIndex))]);
ylim([0,2*max(sparams.totalTime)]);
title('Fidelity: $$|\langle\Psi_0(t_j)|\Psi_{\rm sim}(t_j)\rangle|^2$$','interpreter','latex','fontsize',15);
view(2);
colormap(jet);
colorbar;

if sparams.calculateStarkShift
    for ii = 1:length(sparams.totalTime)
        if mod(ii-1,6) == 0
            fig = figure;
            pos = get(fig,'position');
            set(fig,'position',[pos(1:2)/4 pos(3)*2.0 pos(4)*1.25]);
        end

        subplot(2,3,mod(ii-1,6)+1);
        hold on;
        plot(sparams.tStarkShift(ii,:),(sparams.vShift(ii,:)-sparams.vShift(ii,1))*1E-6);
        plot(sparams.tStarkShift(ii,:),(sparams.vShiftGround(ii,:)-sparams.vShift(ii,1))*1E-6);
        xlabel('Time index [s]','interpreter','latex','fontsize',10);
        ylabel('$\nu - \nu_0$ [MHz]','interpreter','latex','fontsize',10);
        title(['Shuttling Simulation ' num2str(sparams.totalTime(ii)) '[s]'],'interpreter','latex','fontsize',10);
        drawnow;
    end
end
%%
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
    
    
    
    