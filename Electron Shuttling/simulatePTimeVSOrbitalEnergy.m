clear;
shuttleParameterFile;
sparams.adiabaticPulseType = 'effective';

tc = (10:10:100)*1E-6*sparams.ee;
orbitalSpacing = logspace(-2,-5,150)*sparams.ee;

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
sparams.nPulsePoints = 1000;
sparams.dt = 5E-14;

nn = 0; % Iteration tracker
for jj = 1:length(tc)
    dBounds = [0, 1.5E-1;0, 1.5E-1]*sparams.ee;
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
            sparams, adiabThresh, dBounds, 1, effHamiltonianParams );
        
        fprintf(1,'Pulse found with time: %0.3E [s]\n',pulseTimeResults(jj,ii));

    %     [fidTemp, ~, ~, ~, ~, ~] = simulateEffectiveShuttling(sparams, xx, detuningPulse, pulseTimeResults(ii), effHamiltonianParams, NaN, 1);
    %     fidelity(ii) = fidTemp(end);
    end
end

save('TEMP.mat');