clear;
shuttleParameterFile;
sparams.adiabaticPulseType = 'effective';
% sparams.potDir = 'C:\Users\bbuonaco\Documents\GitHub\Simulated Potentials\orbitalSpacingOutputs\';
sparams.potDir = 'C:\Users\bbuonaco\Documents\GitHub\Simulated Potentials\orbitalSpacingOutputs - Correct\';
dotSizes = 30:10:650;

gapSize = 30E-9;
color = {'b','r','g','m','y'};
orbitalSpacing = zeros(1,length(dotSizes));

for ii = 1:length(dotSizes)
    currFolderFile = sprintf('QM_ON_GWIDTHX_%d\\output\\Quantum\\wf_energy_spectrum_quantum_region_X3_0000.dat',dotSizes(ii));
    fid = fopen([sparams.potDir currFolderFile],'r');
    tline = fgetl(fid);
    tline = fgetl(fid);
    energy0 = str2num(tline);
    energy0 = energy0(2);
    tline = fgetl(fid);
    energy1 = str2num(tline);
    energy1 = energy1(2);
    fclose(fid);
    
    orbitalSpacing(ii) = (energy1 - energy0)*sparams.ee;
end

% orbitalInterpolant = griddedInterpolant({dotSizes},orbitalSpacing,'spline');
% dotSizes = linspace(min(dotSizes),max(dotSizes),length(dotSizes));
% orbitalSpacing = orbitalInterpolant({dotSizes});

figure;
hold on;
% plot(dotSizes,orbitalSpacing/sparams.ee,'Linewidth',2);
plot(dotSizes,orbitalSpacing/sparams.ee,'Linewidth',1.5);
%%
tc = linspace(10,100,6)*1E-6*sparams.ee;
% tc = [10,100]*1E-6*sparams.ee;
% tc = 100*1E-6*sparams.ee;
orbitalSpacing = logspace(-2.5,-5,200)*sparams.ee;
% orbitalSpacing = [10^-3,10^-5]*sparams.ee;
% orbitalSpacing = orbitalSpacing(end);

% Get an adiabatic voltage pulse based on the voltage
% points given earlier and using the effective
% Hamiltonian parameters
sparams.includeSpin = 0;
sparams.includeValley = 0;
sparams.includeT2 = 0;
sparams.includeExcitedOrbital = 1;
sparams.adiabaticPulseType = 'effective';

adiabThresh = [0.005,0.005];
pulseTimeResults = zeros(length(tc),length(orbitalSpacing));
fidelity = zeros(length(tc),length(orbitalSpacing));
sparams.nPulsePoints = 600;
sparams.dt = 5E-14;

% tc = 1E-4*sparams.ee;
% orbitalSpacing = 1E-5*sparams.ee;
nn = 0; % Iteration tracker
for jj = 1:length(tc)
    dBounds = [0, 500*tc(jj);0, 500*tc(jj)];
    for ii = 1:length(orbitalSpacing)
        nn = nn + 1;
        fprintf(1,'(%d/%d) Orbital energy = %.3E [eV], tc = %0.3E [eV]\n',...
            nn,length(orbitalSpacing)*length(tc),orbitalSpacing(ii)/sparams.ee,tc(jj)/sparams.ee);

        t1 = tc(jj);
        t2 = tc(jj);
        t3 = tc(jj);
        t4 = tc(jj);    
        
        % Put all the hamiltonian parameters into a single
        % variable for ease
        effHamiltonianParams = buildEffHamiltonianParamVariable(dBounds(1,:), dBounds(2,:),...
            [t1, t2, t3, t4], 0, 0, 0, 0, 0, 0, orbitalSpacing(ii), orbitalSpacing(ii));

        [sparams, detuningPulse, pulseTimeResults(jj,ii)] = getDetuningPulseAdiabatic(...
            sparams, [0.005,0.005], dBounds, 1, effHamiltonianParams );
        
        fprintf(1,'Pulse found with time: %0.3E [s]\n',pulseTimeResults(jj,ii));

    %     [fidTemp, ~, ~, ~, ~, ~] = simulateEffectiveShuttling(sparams, xx, detuningPulse, pulseTimeResults(ii), effHamiltonianParams, NaN, 1);
    %     fidelity(ii) = fidTemp(end);
    end
end