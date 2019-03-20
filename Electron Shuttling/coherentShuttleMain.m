% Load all the parameters for the simulation
clear sparams xx vv pp;
shuttleParameterFile;

% profile on
fprintf(1,'Loading potentials...\n');
[sparams,xx,zz] = loadPotentials(sparams, 1, 'nextnano', 'x', 1);
% profile off
% profile viewer

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
%%
% Check that the potentials were loaded correctly and that the interpolants
% were correctly assembled
debugHere = 0;
if debugHere
    checkPotentialLoad(sparams,xx,zz);
end

testPot = sparams.P2DEGInterpolant({0.2,0.30,0.30,xx});
testPot = squeezeFast(length(sparams.interpableGates),testPot);
plot(xx,testPot/sparams.ee);
%%
% Get our desired votlage pulse.
% vec = [[0.01,0.01];[0.001,0.025];[0.001,0.05]];
vec = [0.05,0.05];
for ii = 1:length(vec(:,1))
    [sparams, vPulse, vPulseTime] = getVoltagePulseAdiabatic(sparams,xx,vec(ii,:),0.8);
end
sparams.voltagePulse = vPulse;

debugHere = 0;
if debugHere
    analyzePulseAdiabicity(sparams,xx,squeeze(vPulse(1,:,:)),...
        vPulseTime(1),1,1:4);
end
%%
%*************************************************************************%
% profile on;

% Using ref https://arxiv.org/pdf/1306.3247.pdf we now find the time
% evolution operator U(t + dT,t) to evolve our initial wavefunction to the
% next wavefunction.  That newly found wavefunction will act as our initial
% state for the next time frame.  This method uses the split operator
% approach
[sparams, sweepVec, ~, ~] = initializeShuttlingSimulation(sparams, pp);

for jj = 1:length(sweepVec)
    if strcmp(sparams.sweptParameter,'time')
        fprintf(1,'Running time sweep shuttling simulation for %.3E (%d/%d)...\n',sweepVec(jj),jj,length(sweepVec));
    elseif strcmp(sparams.sweptParameter,'adiabicity')
        fprintf(1,'Running adiabatic sweep shuttling simulation for [%.3E,%.3E] (%d/%d)...\n',sweepVec(jj,1),sweepVec(jj,2),jj,length(sweepVec));
    end
    
    tic;
    
    fprintf(1,'Building voltage pulse...\n');
    % First step is to get the voltage pulse depending on our simulation
    % type
    if strcmp(sparams.sweptParameter,'time') 
        [sparams, sparams.voltagePulse(jj,:,:)] = getVoltagePulse(sparams,xx);
    elseif strcmp(sparams.sweptParameter,'adiabicity')
        [sparams, sparams.voltagePulse(jj,:,:), sparams.totalTime(jj)] =...
            getVoltagePulseAdiabatic(sparams,xx,sweepVec(jj,:),sparams.voltagePulseBounds,1,NaN);
    end
    
    sparams = simulateCoherentShuttling(sparams, sweepVec, xx, zz, pp, jj, 'split-operator');
    
    % Save the simulation results so far (in case of a shut down mid
    % simulation)
    save([sparams.saveDir sparams.saveFolder 'simulationData'],'sparams','xx','zz','pp');
    
    toc;
end

% profile off
% profile viewer
%*************************************************************************%
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
    
    
    