fprintf(1,'************************\n');
fprintf(1,'* Loading potential... *\n');
fprintf(1,'************************\n\n');

sparams.unitsType = 'Rydberg';
% sparams.unitsType = 'SI - J';
% sparams.unitsType = 'SI - meV';
simparams;

if strcmp(sparams.unitsType,'Rydberg')
    effaB = 1; % [aB*]
    effRy = 1; % [Ry*]
% elseif strcmp(sparams.unitsType,'SI - meV')
%     effaB = 4*pi*sparams.eps*(sparams.hbar^2)/(sparams.me*(sparams.ee^2)); % [m]
%     effRy = (sparams.ee^2)/(8*pi*sparams.eps*effaB); % [meV]
elseif strcmp(sparams.unitsType,'SI - J')
    effaB = 4*pi*sparams.eps*(sparams.hbar^2)/(sparams.me*(sparams.ee^2)); % [m]
    effRy = (sparams.ee^2)/(8*pi*sparams.eps*effaB); % [J]
end

% Scale the dot locations according to our units
sparams.dotLocations = sparams.dotLocations*effaB;

% Make the potential wells
Vimin = 10;
Vimin = Vimin*effRy;
di = 2.3*effaB;
if strcmp(sparams.unitsType,'Rydberg')
    omega = sqrt(Vimin/di^2);
elseif strcmp(sparams.unitsType,'SI - meV') || strcmp(sparams.unitsType,'SI - J')
    omega = sqrt(2*Vimin/(sparams.me*di^2));
end

% Fill in grid parameters
gparams.ngridx = 250;
gparams.ngridy = 250;
gparams.xx = linspace(-8,8,gparams.ngridx)*effaB;
gparams.yy = linspace(-8,8,gparams.ngridy)*effaB;
[gparams.XX,gparams.YY] = meshgrid(gparams.xx,gparams.yy);

gparams.VV = zeros(gparams.ngridy,gparams.ngridx);
for ii = 1:sparams.nDots
    gparams.VV = gparams.VV + -Vimin*exp(-((gparams.XX - sparams.dotLocations(ii,1)).^2 +...
        (gparams.YY - sparams.dotLocations(ii,2)).^2)/di^2);
    sparams.fittedPotentialParameters(ii,:) = [omega, Vimin, sparams.dotLocations(ii,:)];
end

debugHere = 1;
if debugHere == 1
    plotMeshgrid(gparams,gparams.VV);
    set(gca,'visible','off')
    export_fig 'UNSW_3dotPot' -m3
end
fprintf(1,'Done!\n');
%%
bias = linspace(0,0.3,10);
% bias = 0;
origVV = gparams.VV;
energies = zeros(3,length(bias));
sparams.maxOriginHOsX = 14;
sparams.maxOriginHOsY = 14;
sparams.nItinerantOrbitals = 3;
sparams.numElectrons = 3;
sparams.spinSubspaces = [2];
sparams.nOutputtedEnergies = 3;

for ii = 1:length(bias)
    fprintf(1,'****************\n');
    fprintf(1,'Sim ind: %d/%d\n', ii, length(bias));
    fprintf(1,'****************\n');
    gparams.VV = origVV -bias(ii)*exp(-((gparams.XX - sparams.dotLocations(1,1)).^2 +...
        (gparams.YY - sparams.dotLocations(1,2)).^2)/di^2);
    
    optOmegaFlag = 0;
    debugFlag = 0;
    
    [eVecs, ens,~,SEens] = calculateManyBodySpectra_2Bases(sparams, gparams, optOmegaFlag, debugFlag, CMEs_lib);
    energies(:,ii) = diag(ens);
end
%%
cm = viridis(256);
energies = energies - energies(1,:);
figure('Color','white');
hold on;
plot(bias*5.93,energies(1,:)*5.93,'Linewidth',2,'Color',cm(20,:));
plot(bias*5.93,energies(2,:)*5.93,'Linewidth',2,'Color',cm(110,:));
plot(bias*5.93,energies(3,:)*5.93,'Linewidth',2,'Color',cm(245,:));
set(gca,'TickLabelInterpreter','latex','Fontsize',14);
xlabel('$V_1$ [meV]','Interpreter','latex','Fontsize',20);
ylabel('Energy [meV]','Interpreter','latex','Fontsize',20);
ylim([-0.05,0.6]);
xlim([min(bias),max(bias)]*5.93);
% plotMeshgrid(gparams,gparams.VV);
export_fig 'UNSW_manyElecSpectra' -m3

%%
% syms J J23
% H = 1/4*[-J23, 2*J23, 2*J; 2*J23, -J23, 2*J; 2*J, 2*J, -2*J + J23];
% simplify(eig(H))

%%
J = 2/3*(energies(3,:) - energies(2,:));
J23 = 1/3*(energies(2,:) + 2*energies(3,:));

figure;
hold on;
plot(bias,J);
plot(bias,J23);

%%
e1 = zeros(1,length(bias));
e2 = J23 - J;
e3 = J23 + J/2;

figure('Color','white');
hold on;
plot(bias*5.93,e1*5.93,'Linewidth',2,'Color',cm(20,:));
plot(bias*5.93,e2*5.93,'Linewidth',2,'Color',cm(110,:));
plot(bias*5.93,e3*5.93,'Linewidth',2,'Color',cm(245,:));
set(gca,'TickLabelInterpreter','latex','Fontsize',14);
xlabel('$V_1$ [meV]','Interpreter','latex','Fontsize',20);
ylabel('Energy [meV]','Interpreter','latex','Fontsize',20);
ylim([-0.05,0.6]);
xlim([min(bias),max(bias)]*5.93);





