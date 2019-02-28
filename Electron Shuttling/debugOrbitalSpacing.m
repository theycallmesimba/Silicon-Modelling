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



