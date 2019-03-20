load('pTimeVSOrbSpacing.mat');
orbSpacing1 = orbitalSpacing;
tc1 = tc;
pTime1 = pulseTimeResults;

load('optimal-geometry-data.mat');
orbSpacing2 = orbitalSpacing;
tc2 = tc;
pTime2 = pulseTimeResults;

load('optimal-geometry-data2-TEMP.mat');
orbSpacing3 = orbitalSpacing;
tc3 = tc;
pTime3 = pulseTimeResults;

for ii = floor(length(tc2)/2):-1:1
    tc2(2*ii) = [];
    pTime2(2*ii,:) = [];
end

for ii = floor(length(tc3)/2):-1:1
    tc3(2*ii) = [];
    pTime3(2*ii,:) = [];
end
%%
orbSpacingAxisMin = max([min(orbSpacing1), min(orbSpacing2), min(orbSpacing3)]);
orbSpacingAxisMax = min([max(orbSpacing1), max(orbSpacing2), max(orbSpacing3)]);

figure;
subplot(2,2,1);
title('1');
loglog(orbSpacing1,pTime1);
xlim([orbSpacingAxisMin,orbSpacingAxisMax]);
subplot(2,2,2);
title('2');
loglog(orbSpacing2,pTime2);
xlim([orbSpacingAxisMin,orbSpacingAxisMax]);
subplot(2,2,3);
title('3');
loglog(orbSpacing3,pTime3);
xlim([orbSpacingAxisMin,orbSpacingAxisMax]);

%%
sparams.includeSpin = 0;
sparams.includeValley = 0;
sparams.includeT2 = 0;
sparams.includeExcitedOrbital = 1;
sparams.adiabaticPulseType = 'effective';
syms tc epsL epsR EoL EoR deltaL deltaR
effHamiltonianParams = buildEffHamiltonianParamVariable(epsL, epsR,...
    [tc, tc, tc, tc], 0, 0, 0, 0, 0, 0, EoL, EoR);
Heff = constructEffectiveHamiltonian( sparams, effHamiltonianParams)

H = updateEffectiveHamiltonianDetuning(Heff,deltaL,deltaR)

%%
load('orbVSpTime-100tc.mat');
orbVSpTime = griddedInterpolant({fliplr(orbitalSpacing)},fliplr(pulseTimeResults));
figure;
plot(orbitalSpacing,pulseTimeResults);

sparams.potDir = 'C:\Users\bbuonaco\Documents\GitHub\Simulated Potentials\orbitalSpacingOutputs - Correct\';
dotSizes = 30:10:650;
gapSize = 30E-9;
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

orbitalSpacing = smooth(orbitalSpacing)';
pTimes = orbVSpTime(orbitalSpacing);
figure;
plot(orbitalSpacing,pTimes);

fraction = linspace(1,2,length(dotSizes));
for ii = 1:length(dotSizes)
    dotSizes(ii) = dotSizes(ii)/fraction(ii);
end
figure;
velocity = (dotSizes*1E-9 + gapSize)./pTimes;
plot(dotSizes,velocity);



