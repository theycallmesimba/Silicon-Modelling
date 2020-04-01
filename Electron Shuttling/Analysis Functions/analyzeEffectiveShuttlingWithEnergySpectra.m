% load('Ez-40-tc-75-SO-2.0-VLvsdPhi.mat');
shuttleParameterFile;
sparams.nStoreDataFrames = 1000;
sparams.nPulsePoints = 500;
sparams.dt = 5E-14;

dPhi =  pi/2;
valleyL = 100*1E-6*sparams.ee;
valleyR = [150]*1E-6*sparams.ee;
spinOrbit = 2.0*1E-6*sparams.ee;
Ez = [40]*1E-6*sparams.ee;
tc = 75E-6*sparams.ee;
EL = 0;
ER = 0;
adiabThresh = 0.005;
dBounds = [0, 1500E-6;0, 1500E-6]*sparams.ee;

sparams.includeExcitedOrbital = 0;
sparams.includeSecondSpin = 1;
sparams.stateIndices = 1;
sparams.nnIndices = 1;

pulseTime = [];
fidelity = zeros(1,sparams.nStoreDataFrames);

effHamiltonianParams = buildEffHamiltonianParamVariable(NaN, NaN,...
    tc, Ez, 0, valleyL*exp(1i*dPhi), valleyR, spinOrbit, spinOrbit, EL, ER);

fprintf(1,'Finding adiabatic pulse...\n');
[sparams, detuningPulse, pTimeTemp] = getDetuningPulseAdiabatic(...
    sparams, [adiabThresh, adiabThresh], dBounds, 1, effHamiltonianParams);

fprintf(1,'Found pulse time = %.6E [s]\n\n',pTimeTemp);
fprintf(1,'Simulating electron shuttling (effective Hamiltonian)...\n');
[rhos, Hams] = simulateEffectiveShuttling(sparams, [], detuningPulse, pTimeTemp, effHamiltonianParams, NaN, 1);