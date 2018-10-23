shuttleParameterFile;
sparams.potDir = 'C:\Users\bbuonaco\Documents\GitHub\Simulated Potentials\Increasing gate size 40nm-500nm\';

gapSize = 20E-9;
dotSizes = linspace(40,500,93);
dotSizes = [dotSizes(1:getClosestArrayIndex(325,dotSizes)-1),dotSizes(getClosestArrayIndex(325,dotSizes)+1:end)];
% dotSizes = [320,325,330];
color = {'b','r','g'};
orbitalSpacing = zeros(1,length(dotSizes));
orbitalSpacingTheor = zeros(1,length(dotSizes));

% figure;
% hold on;
for ii = 1:length(dotSizes)
    currFile = sprintf('GATE_SIZE_%d.csv',dotSizes(ii));
    [xx, zz, currPot2D] = loadPotentialFile([sparams.potDir, currFile]);
    xx = xx*1E-9;
    zz = zz*1E-9;
    currPot2D = -sparams.ee*currPot2D;
    
    currPot2DEG = currPot2D(getClosestArrayIndex(-0.5*1E-9, zz),:);
    [~, ens] = solve1DSingleElectronSE(sparams, 2, xx, currPot2DEG);

    orbitalSpacing(ii) = ens(2,2) - ens(1,1);
    
%     plot(xx,currPot2DEG/sparams.ee,'Color',color{ii});
%     line([min(xx),max(xx)],[ens(1,1),ens(1,1)]/sparams.ee,'Color',color{ii});
%     line([min(xx),max(xx)],[ens(2,2),ens(2,2)]/sparams.ee,'Color',color{ii});
%     line([min(xx),max(xx)],[ens(3,3),ens(3,3)]/sparams.ee,'Color',color{ii});
    
    orbitalSpacingTheor(ii) = 3*pi^2*sparams.hbar^2/(2*sparams.me*(dotSizes(ii)*1E-9)^2);
end

% orbitalInterpolant = griddedInterpolant({dotSizes},orbitalSpacing,'spline');
% dotSizes = linspace(40,500,200);
% orbitals = orbitalInterpolant({dotSizes});
orbitals = orbitalSpacing;
% figure;
% hold on;
% plot(dotSizes,orbitals/sparams.ee);
% % plot(linspace(40,500,93),orbitalSpacingTheor)

tc = 50E-6*sparams.ee;
t1 = tc;
t2 = tc;
t3 = tc;
t4 = tc;
% Get an adiabatic voltage pulse based on the voltage
% points given earlier and using the effective
% Hamiltonian parameters
sparams.includeSpin = 0;
sparams.includeValley = 0;
sparams.includeT2 = 0;
sparams.includeExcitedOrbital = 1;

adiabThresh = [0.005,0.005];
dBounds = [0, 5000E-6;0, 5000E-6]*sparams.ee;
pulseTimeResults = zeros(1,length(orbitals));
fidelity = zeros(1,length(orbitals));
sparams.nPulsePoints = 500;
sparams.dt = 5E-14;

for ii = 1:length(orbitals)
% for ii = 200
    fprintf(1,'(%d/%d) Orbital energy = %.3E [eV], tc = %0.3E [eV]\n',...
        ii,length(orbitals),orbitals(ii)/sparams.ee,tc/sparams.ee);

    % Put all the hamiltonian parameters into a single
    % variable for ease
    effHamiltonianParams = buildEffHamiltonianParamVariable(dBounds(1,:), dBounds(2,:),...
        [t1, t2, t3, t4], 0, 0, 0, 0, 0, 0, orbitals(ii), orbitals(ii));

    [sparams, detuningPulse, pulseTimeResults(ii)] = getDetuningPulseAdiabatic(...
        sparams, [0.005,0.005], dBounds, 1, effHamiltonianParams );
    
    [fidTemp, ~, ~, ~, ~, ~] = simulateEffectiveShuttling(sparams, xx, detuningPulse, pulseTimeResults(ii), effHamiltonianParams, NaN, 1);
    fidelity(ii) = fidTemp(end);
end
%%
velocity = (dotSizes*1E-9 + gapSize)./pulseTimeResults;
figure;
yyaxis left
plot(dotSizes,velocity);
yyaxis right
plot(dotSizes,fidelity);

%%
ii = 27;
sparams.dt = 5E-14;
tc = 50E-6*sparams.ee;
t1 = tc;
t2 = tc;
t3 = tc;
t4 = tc;
effHamiltonianParams = buildEffHamiltonianParamVariable(dBounds(1,:), dBounds(2,:),...
    [t1, t2, t3, t4], 0, 0, 0, 0, 0, 0, orbitals(ii), orbitals(ii));

[sparams, detuningPulse, pT] = getDetuningPulseAdiabatic(...
    sparams, [0.05,0.05], dBounds, 1, effHamiltonianParams );

[f, ~, ~, ~, ~, ~] = simulateEffectiveShuttling(sparams, xx, detuningPulse, pT, effHamiltonianParams, NaN, 1);
fprintf(1,'Pulse Time %.3E with fidelity %0.6E\n',pT,f(end));


