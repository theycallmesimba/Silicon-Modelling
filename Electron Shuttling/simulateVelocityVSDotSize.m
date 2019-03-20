shuttleParameterFile;
% sparams.potDir = 'C:\Users\bbuonaco\Documents\GitHub\Simulated Potentials\Increasing gate size 40nm-500nm\';
% dotSizes = linspace(40,500,93);
% dotSizes = [dotSizes(1:getClosestArrayIndex(325,dotSizes)-1),dotSizes(getClosestArrayIndex(325,dotSizes)+1:end)];
% dotSizes = [40,60,80];
% sparams.potDir = 'C:\Users\bbuonaco\Documents\GitHub\Simulated Potentials\Increasing gate size 320-330nm & 150-190nm\';
% dotSizes = linspace(320,330,11);
% dotSizes = linspace(150,190,41);
% sparams.potDir = 'C:\Users\bbuonaco\Documents\GitHub\Simulated Potentials\New Grid settings Increasing gate size 150-190\';
% dotSizes = linspace(150,190,21);
sparams.potDir = 'C:\Users\bbuonaco\Documents\GitHub\Simulated Potentials\Expanding dot size- Fine Grid\';
dotSizes = linspace(40,400,181);

gapSize = 20E-9;
color = {'b','r','g','m','y'};
orbitalSpacing = zeros(1,length(dotSizes));
orbitalSpacingTheor = zeros(1,length(dotSizes));

debugFlag = 0;
% dotSizes = 150;
for ii = 1:length(dotSizes)
%     currFile = sprintf('GATE_SIZE_%d.csv',dotSizes(ii));
    currFile = sprintf('GATE_WIDTH_X_%d.csv',dotSizes(ii));
    [xx, zz, currPot2D] = loadPotentialFile(sparams, [sparams.potDir, currFile],...
        1,'csv','xz',0);
    
    currPot2DEG = currPot2D(getClosestArrayIndex(-0.5*1E-9, zz),:);
    [wfs, ens] = solve1DSingleElectronSE(sparams, 3, xx, currPot2DEG);

    orbitalSpacing(ii) = ens(2,2) - ens(1,1);
    
    if debugFlag
        if ii == 1
            figure;
            hold on;
        end
        normConst = 2E7*100;
        plot(xx/1E-9,(wfs(:,1)).^2/normConst);
        plot(xx/1E-9,(wfs(:,2)).^2/normConst);
        plot(xx/1E-9,(currPot2DEG - min(currPot2DEG))/sparams.ee,'Color',color{ii},'linewidth',2);
%         line([min(xx),max(xx)]/1E-9,[ens(1,1),ens(1,1)]/sparams.ee,'Color',color{ii},'linewidth',2);
%         line([min(xx),max(xx)]/1E-9,[ens(2,2),ens(2,2)]/sparams.ee,'Color',color{ii},'linewidth',2);
    %     line([min(xx),max(xx)],[ens(3,3),ens(3,3)]/sparams.ee,'Color',color{ii});
    end
    
    orbitalSpacingTheor(ii) = 3*pi^2*sparams.hbar^2/(2*sparams.me*(dotSizes(ii)*1E-9)^2);
end
if debugFlag
    set(gca,'Fontsize',14);
    xlabel('Position [nm]','Fontsize',18,'Interpreter','latex');
    ylabel('Energy [eV]','Fontsize',18,'Interpreter','latex');
end

% orbitalInterpolant = griddedInterpolant({dotSizes},orbitalSpacing,'spline');
% dotSizes = linspace(40,500,200);
% orbitals = orbitalInterpolant({dotSizes});
% orbitals = orbitalSpacing;
figure;
hold on;
plot(dotSizes,orbitalSpacing/sparams.ee,'Linewidth',2);
% plot(dotSizes,orbitalSpacingTheor/sparams.ee);
%%
% tc = linspace(10,100,11)*1E-6*sparams.ee;
tc = linspace(10,100,6)*1E-6*sparams.ee;
orbitalSpacing = [4.84289E-22,6.79833E-24];

% Get an adiabatic voltage pulse based on the voltage
% points given earlier and using the effective
% Hamiltonian parameters
sparams.includeSpin = 0;
sparams.includeValley = 0;
sparams.includeT2 = 0;
sparams.includeExcitedOrbital = 1;
sparams.adiabaticPulseType = 'effective';

adiabThresh = [0.005,0.005];
% dBounds = [0, 4500E-6;0, 4500E-6]*sparams.ee;
dBounds = [0, 5000E-6;0, 5000E-6]*sparams.ee;
pulseTimeResults = zeros(length(tc),length(orbitalSpacing));
fidelity = zeros(length(tc),length(orbitalSpacing));
sparams.nPulsePoints = 500;
sparams.dt = 5E-14;
nn = 0; % Iteration tracker
for jj = 1:length(tc)
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
%%
shuttleParameterFile;

sparams.dt = 5E-14;
% tc = 31.3E-6*sparams.ee;
valleyL = valleyL(11);
valleyR = [100]*1E-6*sparams.ee;
spinOrbit = spinOrbit(3);
Ez = [40]*1E-6*sparams.ee;
Ex = [0]*1E-6*sparams.ee;
dBounds = [0, 1500E-6;0, 1500E-6]*sparams.ee;







