% Zeroth step, load the potentials
shuttleParameterFile;

fprintf(1,'Loading potentials...\n');
[sparams,xx,zz] = loadPotentials(sparams);
sparams.nxGrid = length(xx);
sparams.nzGrid = length(zz);
sparams.dx = xx(2) - xx(1);
sparams.dz = zz(2) - zz(1);

% Find which index corresponds to where the 2DEG should be
sparams.twoDEGindZ = getClosestArrayIndex(-0.5*1E-9, zz);
for ii = 1:length(sparams.potentials)
    sparams.potentials(ii).pot2DEG = sparams.potentials(ii).pot2D(sparams.twoDEGindZ,:);
end

% Now we want to make the potential interpolant object (both 2D and 2DEG)
sparams.interpType = 'spline';
sparams.extrapType = 'spline';
sparams = makePotentialsInterpolants(sparams,xx,zz);

% First step is to define our voltage search range and then find how the
% detuning changes in that window.
adiabThresh = [0.005,0.005];
vBounds = [0.54,0.6; 0.54,0.6; 0.4,0.6];

[~, voltagePulse, ~] = getVoltagePulseAdiabatic( sparams, xx, adiabThresh, vBounds, 0, NaN );
[epsL, epsR] = getDetuningVsVoltagePulse( sparams, xx, voltagePulse, 1 );

% Second step is to get the tunnel coupling at the maximum voltage point in
% our pulse
tcVArg = voltagePulse(:,1);
tcVArg(sparams.gatesUsedInPulse(2)) = max(voltagePulse(sparams.gatesUsedInPulse(2),:));
currPotential = sparams.P2DEGInterpolant(getInterpolantArgument(tcVArg,xx));
currPotential = squeezeFast(sparams.numOfGates,currPotential)';
tc = calculateTunnelCoupling( sparams, xx, currPotential );
%%
%*********************************************************************%
%%
%*COMPLETE*
saveName = 'Ez-75-tc-40-dPhi-2pi5-VLvsSO';
dPhi =  2*pi/5;
valleyL = linspace(25,250,50)*1E-6*sparams.ee;
valleyR = [150]*1E-6*sparams.ee;
spinOrbit = linspace(0,2.0,50)*1E-6*sparams.ee;
Ex = [0]*1E-6*sparams.ee;
Ez = [75]*1E-6*sparams.ee;
tc = 40E-6*sparams.ee;
EL = 0;
ER = 0;
%%
% COMPLETE
saveName = 'Ez-75-tc-40-dPhi-3pi5-VLvsSO';
dPhi =  3*pi/5;
valleyL = linspace(25,250,50)*1E-6*sparams.ee;
valleyR = [150]*1E-6*sparams.ee;
spinOrbit = linspace(0,2.0,50)*1E-6*sparams.ee;
Ex = [0]*1E-6*sparams.ee;
Ez = [75]*1E-6*sparams.ee;
tc = 40E-6*sparams.ee;
EL = 0;
ER = 0;
%%
% COMPLETE
saveName = 'Ez-40-tc-75-dPhi-2pi5-VLvsSO';
dPhi =  2*pi/5;
valleyL = linspace(25,250,50)*1E-6*sparams.ee;
valleyR = [150]*1E-6*sparams.ee;
spinOrbit = linspace(0,2.0,50)*1E-6*sparams.ee;
Ex = [0]*1E-6*sparams.ee;
Ez = [40]*1E-6*sparams.ee;
tc = 75E-6*sparams.ee;
EL = 0;
ER = 0;
%%
% COMPLETE
saveName = 'Ez-40-tc-75-dPhi-3pi5-VLvsSO';
dPhi =  3*pi/5;
valleyL = linspace(25,250,50)*1E-6*sparams.ee;
valleyR = [150]*1E-6*sparams.ee;
spinOrbit = linspace(0,2.0,50)*1E-6*sparams.ee;
Ex = [0]*1E-6*sparams.ee;
Ez = [40]*1E-6*sparams.ee;
tc = 75E-6*sparams.ee;
EL = 0;
ER = 0;
%%
% COMPLETE
saveName = 'Ez-75-tc-40-SO-0.4-VLvsdPhi.mat';
dPhi =  linspace(0,pi,50);
valleyL = linspace(25,250,50)*1E-6*sparams.ee;
valleyR = [150]*1E-6*sparams.ee;
spinOrbit = 0.4*1E-6*sparams.ee;
Ex = [0]*1E-6*sparams.ee;
Ez = [75]*1E-6*sparams.ee;
tc = 40E-6*sparams.ee;
EL = 0;
ER = 0;
%%
% COMPLETE
saveName = 'Ez-75-tc-40-SO-2.0-VLvsdPhi.mat';
dPhi =  linspace(0,pi,50);
valleyL = linspace(25,250,50)*1E-6*sparams.ee;
valleyR = [150]*1E-6*sparams.ee;
spinOrbit = 2.0*1E-6*sparams.ee;
Ex = [0]*1E-6*sparams.ee;
Ez = [75]*1E-6*sparams.ee;
tc = 40E-6*sparams.ee;
EL = 0;
ER = 0;
%%
% COMPLETE
saveName = 'Ez-40-tc-75-SO-0.4-VLvsdPhi.mat';
dPhi =  linspace(0,pi,50);
valleyL = linspace(25,250,50)*1E-6*sparams.ee;
valleyR = [150]*1E-6*sparams.ee;
spinOrbit = 0.4*1E-6*sparams.ee;
Ex = [0]*1E-6*sparams.ee;
Ez = [40]*1E-6*sparams.ee;
tc = 75E-6*sparams.ee;
EL = 0;
ER = 0;
%%
% RUNNING
saveName = 'Ez-40-tc-75-SO-2.0-VLvsdPhi.mat';
dPhi =  linspace(0,pi,50);
valleyL = linspace(25,250,50)*1E-6*sparams.ee;
valleyR = [150]*1E-6*sparams.ee;
spinOrbit = 2.0*1E-6*sparams.ee;
Ex = [0]*1E-6*sparams.ee;
Ez = [40]*1E-6*sparams.ee;
tc = 75E-6*sparams.ee;
EL = 0;
ER = 0;
%%
% saveName = 'TEMP';
T2Sweep = 100000E-9;
adiabThresh = 0.005;
vBounds = [0.44,0.60; 0.44,0.6; 0.44,0.6];
dBounds = [0, 1500E-6;0, 1500E-6]*sparams.ee;
sparams.dt = 5E-14;

sparams.includeExcitedOrbital = 0;
sparams.includeSecondSpin = 1;
sparams.stateIndices = 1;
sparams.nnIndices = 1;

pulseTime = zeros(length(adiabThresh),length(T2Sweep),length(valleyL),length(valleyR),...
    length(dPhi),length(spinOrbit),length(Ez),length(Ex));
fidelity = zeros(length(adiabThresh),length(T2Sweep),length(valleyL),length(valleyR),...
    length(dPhi),length(spinOrbit),length(Ez),length(Ex),sparams.nStoreDataFrames);

totPurity = zeros(length(adiabThresh),length(T2Sweep),length(valleyL),length(valleyR),...
    length(dPhi),length(spinOrbit),length(Ez),length(Ex),sparams.nStoreDataFrames);
orbPurity = zeros(length(adiabThresh),length(T2Sweep),length(valleyL),length(valleyR),...
    length(dPhi),length(spinOrbit),length(Ez),length(Ex),sparams.nStoreDataFrames);
valPurity = zeros(length(adiabThresh),length(T2Sweep),length(valleyL),length(valleyR),...
    length(dPhi),length(spinOrbit),length(Ez),length(Ex),sparams.nStoreDataFrames);
spinPurity = zeros(length(adiabThresh),length(T2Sweep),length(valleyL),length(valleyR),...
    length(dPhi),length(spinOrbit),length(Ez),length(Ex),sparams.nStoreDataFrames);

orbExpectation = zeros(length(adiabThresh),length(T2Sweep),length(valleyL),length(valleyR),...
    length(dPhi),length(spinOrbit),length(Ez),length(Ex),3,sparams.nStoreDataFrames);
valExpectation = zeros(length(adiabThresh),length(T2Sweep),length(valleyL),length(valleyR),...
    length(dPhi),length(spinOrbit),length(Ez),length(Ex),3,sparams.nStoreDataFrames);
spinExpectation = zeros(length(adiabThresh),length(T2Sweep),length(valleyL),length(valleyR),...
    length(dPhi),length(spinOrbit),length(Ez),length(Ex),3,sparams.nStoreDataFrames);


uu = 1;
mm = 0;
startTime = clock;
nLoops = length(T2Sweep)*length(valleyL)*length(valleyR)*length(dPhi)*...
    length(spinOrbit)*length(Ez)*length(Ex)*length(adiabThresh);
for cc = 1:length(adiabThresh)
    for ii = 1:length(T2Sweep)
        for jj = 1:length(valleyL)
        for kk = 1:length(valleyR)
        for pp = 1:length(dPhi)
            % If we don't need valley, just omit from the simulation as it
            % messes with the adiabatic pulse finding due to the degenerate
            % ground states
            if valleyL(jj) == 0 || valleyR(kk) == 0
                sparams.includeValley = 0;
            else
                sparams.includeValley = 1;
            end

            for aa = 1:length(spinOrbit)
                for bb = 1:length(Ez)
                for dd = 1:length(Ex)
                    % If we don't need spin, just omit it from the
                    % simulation as it messes with the adiabatic pulse
                    % finding due to the degenerate ground states
                    if Ez(bb) == 0 && Ex(dd) == 0
                        sparams.includeSpin = 0;
                    else
                        sparams.includeSpin = 1;
                    end

                    if mm ~= 0
                        % Estaimte the remaining runtime
                        elapsedTime = etime(clock,startTime);
                        
                        avgTimePerLoop = elapsedTime/mm;
                        avgTimeRemaining = avgTimePerLoop*(nLoops - mm);
                        if avgTimeRemaining > 3600
                            fprintf(1,'\nTime remaining: %4.1f hours\n\n',avgTimeRemaining/3600);
                        elseif avgTimeRemaining > 60
                            fprintf(1,'\nTime remaining: %4.1f minutes\n\n',avgTimeRemaining/60);
                        else
                            fprintf(1,'\nTime remaining: %4.1f seconds\n\n',avgTimeRemaining);
                        end
                    end
                    
                    % Display to use what parameters are bein probed right
                    % now
                    mm = mm + 1;
                    fprintf(1,'********************************************\n');
                    fprintf(1,'(%d/%d): Adiabatic Threshold = %.3E [arb]\n',...
                        mm, nLoops, adiabThresh(cc));
                    fprintf(1,'tc = %.3E [eV], T2 = %.3E [s]\n',tc/sparams.ee,T2Sweep(ii));
                    fprintf(1,'ValleyL = %.3E [eV], ValleyR = %.3E [eV], Phase = %.3f\n',...
                        valleyL(jj)/sparams.ee, valleyR(kk)/sparams.ee, dPhi(pp));
                    fprintf(1,'Ez = %.3E [eV], Ex = %.3E [eV], Spin Orbit = %.3E [eV]\n\n',...
                        Ez(bb)/sparams.ee, Ex(dd)/sparams.ee, spinOrbit(aa)/sparams.ee);

                    % Put all the hamiltonian parameters into a single
                    % variable for ease
                    effHamiltonianParams = buildEffHamiltonianParamVariable(epsL, epsR,...
                        tc, Ez(bb), Ex(dd), valleyL(jj)*exp(1i*dPhi(pp)), valleyR(kk), spinOrbit(aa), spinOrbit(aa), EL, ER);

                    % Get an adiabatic voltage pulse based on the voltage
                    % points given earlier and using the effective
                    % Hamiltonian parameters
%                     [sparams, voltagePulse, pTimeTemp] = getVoltagePulseAdiabatic(...
%                         sparams, xx, [adiabThresh(cc),adiabThresh(cc)], vBounds, 1, effHamiltonianParams );
%                     fprintf(1,'Pulse found!\nPulse time = %.6E [s]\n',pTimeTemp);
% 
%                     % Using the adiabatic voltage pulse, run an effective
%                     % shuttling simulation.
%                     fprintf(1,'Now doing effective shuttling simulation using found pulse...\n\n');
%                     [rhos, Hams] = simulateEffectiveShuttling(sparams, xx,...
%                         voltagePulse, pTimeTemp, effHamiltonianParams, T2Sweep(ii), 0);
                    
                    %****************%
                    fprintf(1,'Finding adiabatic pulse...\n');
                    [sparams, detuningPulse, pTimeTemp] = getDetuningPulseAdiabatic(...
                        sparams, [adiabThresh(cc), adiabThresh(cc)], dBounds, 1, effHamiltonianParams );
                    fprintf(1,'Found pulse time = %.6E [s]\n\n',pTimeTemp);
                    fprintf(1,'Simulating electron shuttling (effective Hamiltonian)...\n');
                    [rhos, Hams] = simulateEffectiveShuttling(sparams, xx, detuningPulse, pTimeTemp, effHamiltonianParams, NaN, 1);
                    fprintf(1,'********************************************\n');
                    %****************%
                    
                    % Analyze the data from the simulation
                    results = analyzeEffectiveShuttlingResults(sparams, rhos, Hams);
                    
                    pulseTime(cc,ii,jj,kk,pp,aa,bb,dd) = pTimeTemp;
                    fidelity(cc,ii,jj,kk,pp,aa,bb,dd,:) = results.fidelity;
                    
                    totPurity(cc,ii,jj,kk,pp,aa,bb,dd,:) = results.totPur;
                    orbPurity(cc,ii,jj,kk,pp,aa,bb,dd,:) = results.orbPur;
                    valPurity(cc,ii,jj,kk,pp,aa,bb,dd,:) = results.valPur;
                    spinPurity(cc,ii,jj,kk,pp,aa,bb,dd,:) = results.spinPur;
                    
                    orbExpectation(cc,ii,jj,kk,pp,aa,bb,dd,:,:) = results.orbExp;
                    valExpectation(cc,ii,jj,kk,pp,aa,bb,dd,:,:) = results.valExp;
                    spinExpectation(cc,ii,jj,kk,pp,aa,bb,dd,:,:) = results.spinExp;   
                end
                end
            end
        end
        end 
        end
    end
end

save(saveName);
fprintf(1,'Total simulation run time: %4.1f seconds\n\n',etime(clock,startTime));
fprintf(1,'********************************************\n');
% profile viewer
% profile off
%%
phaseInd = 1;
velocity = (40E-9+20E-9)./pulseTime;
plotEffHamResults2D(valleyL, 'valleyL', spinOrbit, 'spinOrbit', pulseTime(1,1,:,1,phaseInd,:,1,1), 'pulseTime');
plotEffHamResults2D(valleyL, 'valleyL', spinOrbit, 'spinOrbit', velocity(1,1,:,1,phaseInd,:,1,1), 'velocity');
plotEffHamResults2D(valleyL, 'valleyL', spinOrbit, 'spinOrbit', fidelity(1,1,:,1,phaseInd,:,1,1,end), 'singlet fidelity');
%%
soInd = 6;
velocity = (40E-9+20E-9)./pulseTime;
cutOffdPhi = 1;
plotEffHamResults2D(valleyL, 'valleyL', dPhi(1:end-cutOffdPhi), 'dPhi', pulseTime(1,1,:,1,1:end-cutOffdPhi,soInd,1,1), 'pulseTime');
plotEffHamResults2D(valleyL, 'valleyL', dPhi(1:end-cutOffdPhi), 'dPhi', velocity(1,1,:,1,1:end-cutOffdPhi,soInd,1,1), 'velocity');
plotEffHamResults2D(valleyL, 'valleyL', dPhi(1:end-cutOffdPhi), 'dPhi', fidelity(1,1,:,1,1:end-cutOffdPhi,soInd,1,1,end), 'singlet fidelity');
% plotEffHamResults2D(valleyL, 'valleyL', dPhi(1:end-cutOffdPhi), 'dPhi', orbExpectation(1,1,:,1,1:end-cutOffdPhi,soInd,1,1,3,end), 'orbital expectation');
% export_fig 'test.jpg' -m5
%%
dP = 0;
vL = [158.3]*1E-6*sparams.ee;
vR = [150]*1E-6*sparams.ee;
sO = [2]*1E-6*sparams.ee;
Ez = [75]*1E-6*sparams.ee;
Ex = [0]*1E-6*sparams.ee;
tc = 40E-6*sparams.ee;
aThresh = 0.005;
sparams.includeExcitedOrbital = 0;
sparams.includeSecondSpin = 1;
sparams.stateIndices = 1;
sparams.nnIndices = 1;
det = 1500E-6;
dBounds = [0, det;0, det]*sparams.ee;

effHamiltonianParams = buildEffHamiltonianParamVariable(epsL, epsR,...
    tc, Ez, Ex, vL*exp(1i*dP), vR, sO, sO, EL, ER);
%****************%
% [sparams, detuningPulse, pTimeTemp] = getDetuningPulseAdiabatic(...
%     sparams, [aThresh,aThresh], dBounds, 1, effHamiltonianParams );
% pTimeTemp
% fprintf(1,'Pulse found!\nPulse time = %.6E [s]\n',pTimeTemp);
% [rhos, Hams] = simulateEffectiveShuttling(sparams, xx, detuningPulse, pTimeTemp, effHamiltonianParams, NaN, 1);
%****************%
% results = analyzeEffectiveShuttlingResults(sparams, rhos, Hams);
% results.fidelity(end)
%%
%%nPts = 400;
dP = pi/3;
vL = [150]*1E-6*sparams.ee;
vR = [150]*1E-6*sparams.ee;
sO = [2]*1E-6*sparams.ee;
Ez = [40]*1E-6*sparams.ee;
Ex = [0]*1E-6*sparams.ee;
tc = 75*1E-6*sparams.ee;
E0 = 0E-3*sparams.ee;
nPts = 500;

sparams.includeSpin = 0;
sparams.includeValley = 1;
sparams.includeSecondSpin = 0;
sparmas.includeExcitedOrbital = 0;

% det = 1500E-6/2;
det = 800E-6/2;
epsL = linspace(-det,det,nPts)*sparams.ee + E0;
epsR = linspace(det,-det,nPts)*sparams.ee + E0;


effHamiltonianParams = buildEffHamiltonianParamVariable(epsL, epsR,...
    tc, Ez, Ex, vL*exp(-1i*dP/2), vR*exp(1i*dP/2), sO, sO, NaN, NaN);

% analyzeEffectiveEnergySpectrum(sparams, effHamiltonianParams, 'detuning', {'spin','Z'});
analyzeEffectiveEnergySpectrum(sparams, effHamiltonianParams, 'detuning', {'none',''});

export_fig 'valley-orbitEnergySpectrum' -m5


