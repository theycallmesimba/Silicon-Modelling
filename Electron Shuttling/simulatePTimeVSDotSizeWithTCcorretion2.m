clear;
shuttleParameterFile;
sparams.potDir = 'C:\Users\bbuonaco\Documents\GitHub\Simulated Potentials\orbitalSpacingOutputs - Correct\';
dotSizesOrb = 30:10:150;

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

% figure;
% plot(dotSizesOrb,orbitalSpacing/sparams.ee);

% Now tc versus dot size
sparams.potDir = 'C:\Users\bbuonaco\Documents\GitHub\Simulated Potentials\tcOutputsVaryingDotandGapSize\';
dotSizesTc = 30:10:150;
gapSizes = 5:5:40;
tc = zeros(length(gapSizes),length(dotSizesTc));
for jj = 1:length(gapSizes)
for ii = 1:length(dotSizesTc)
    currFolderFile = sprintf('QM_ON_GWIDTHX_%d_GAP_%d\\output\\Quantum\\wf_energy_spectrum_quantum_region_X3_0000.dat',dotSizesTc(ii),gapSizes(jj));
    fid = fopen([sparams.potDir currFolderFile],'r');
    tline = fgetl(fid);
    tline = fgetl(fid);
    energy0 = str2num(tline);
    energy0 = energy0(2);
    tline = fgetl(fid);
    energy1 = str2num(tline);
    energy1 = energy1(2);
    fclose(fid);
    
    tempTc = (energy1 - energy0)/2*sparams.ee;
    % Cutoff results with too large of a calculated tc since those are most
    % likely not well defined tcs
    if tempTc >= 2E-4*sparams.ee
        tempTc = NaN;
    end
    tc(jj,ii) = tempTc;
end
end

figure('Color','white');
[XX,YY] = meshgrid(dotSizesTc,gapSizes);
h = surf(XX,YY,tc/sparams.ee);
set(h,'EdgeColor','none');
view(2);
% set(gca,'colorScale','log','Fontsize',16,'TickLabelInterpreter','latex');
set(gca,'colorScale','log','Fontsize',18,'TickLabelInterpreter','latex');
xlabel('$D$ [nm]','FontSize',20,'Interpreter','latex');
ylabel('$G$ [nm]','FontSize',20,'Interpreter','latex');
% shading(gca,'interp');
cb = colorbar;
cb.Ticks = [1E-7,1E-6,1E-5,1E-4];
cb.TickLabelInterpreter = 'latex';
ylabel(cb,'$t_c$ [eV]','FontSize',20,'Interpreter','latex');
pos = get(gca,'Position');
pos(3) = 0.98*pos(3);
set(gca,'Position',pos);
%%
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

