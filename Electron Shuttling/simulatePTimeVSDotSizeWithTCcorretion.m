clear;
shuttleParameterFile;
sparams.potDir = 'C:\Users\bbuonaco\Documents\GitHub\Simulated Potentials\orbitalSpacingOutputs - Correct\';
dotSizesOrb = 30:10:150;
gapSize = 30E-9;

% First get orbital spacing versus dot size
orbitalSpacing = zeros(1,length(dotSizesOrb));
for ii = 1:length(dotSizesOrb)
    currFolderFile = sprintf('QM_ON_GWIDTHX_%d\\output\\Quantum\\wf_energy_spectrum_quantum_region_X3_0000.dat',dotSizesOrb(ii));
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
orbitalSpacing = smooth(orbitalSpacing);

figure;
plot(dotSizesOrb,orbitalSpacing/sparams.ee);

% Now tc versus dot size
sparams.potDir = 'C:\Users\bbuonaco\Documents\GitHub\Simulated Potentials\tunnnelCouplingOutputs - Correct\';
dotSizesTc = 30:10:150;
tc = zeros(1,length(dotSizesTc));
for ii = 1:length(dotSizesTc)
    currFolderFile = sprintf('QM_ON_GWIDTHX_%d\\output\\Quantum\\wf_energy_spectrum_quantum_region_X3_0000.dat',dotSizesTc(ii));
    fid = fopen([sparams.potDir currFolderFile],'r');
    tline = fgetl(fid);
    tline = fgetl(fid);
    energy0 = str2num(tline);
    energy0 = energy0(2);
    tline = fgetl(fid);
    energy1 = str2num(tline);
    energy1 = energy1(2);
    fclose(fid);
    
    tc(ii) = (energy1 - energy0)/2*sparams.ee;
end
tc = smooth(tc);

figure;
plot(dotSizesTc,tc/sparams.ee);

% Get an adiabatic voltage pulse based on the voltage
% points given earlier and using the effective
% Hamiltonian parameters
sparams.includeSpin = 0;
sparams.includeValley = 0;
sparams.includeT2 = 0;
sparams.includeExcitedOrbital = 1;
sparams.adiabaticPulseType = 'effective';

adiabThresh = [0.005,0.005];
pulseTimeResults = zeros(1,length(dotSizesOrb));
sparams.nPulsePoints = 750;
sparams.dt = 5E-14;

nn = 0; % Iteration tracker
for ii = 1:length(dotSizesOrb)
    tcCurr = tc(ii);
    orbSpacCurr = orbitalSpacing(ii);
    dBounds = [0, 750*tcCurr;0, 750*tcCurr];
    
    nn = nn + 1;
    fprintf(1,'(%d/%d) Orbital energy = %.3E [eV], tc = %0.3E [eV]\n',...
        nn,length(dotSizesOrb),orbSpacCurr/sparams.ee,tcCurr/sparams.ee);

    t1 = tcCurr;
    t2 = tcCurr;
    t3 = tcCurr;
    t4 = tcCurr;    

    % Put all the hamiltonian parameters into a single
    % variable for ease
    effHamiltonianParams = buildEffHamiltonianParamVariable(dBounds(1,:), dBounds(2,:),...
        [t1, t2, t3, t4], 0, 0, 0, 0, 0, 0, orbSpacCurr, orbSpacCurr);

    [sparams, ~, pulseTimeResults(ii)] = getDetuningPulseAdiabatic(...
        sparams, adiabThresh, dBounds, 1, effHamiltonianParams );

    fprintf(1,'Pulse found with time: %0.3E [s]\n',pulseTimeResults(ii));
end

velocity = (dotSizesOrb*1E-9 + gapSize)./pulseTimeResults;
figure;
plot(dotSizesOrb,velocity);

