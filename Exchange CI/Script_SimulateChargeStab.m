%% Load nextnano simulations
% Zeroth step, load the potentials
shuttleParameterFile;

aBstar = 2.16; % Convert abStar to nm
Rystar = 42.7; % Convert Ry* to meV

fprintf(1,'Loading potentials...\n');
[sparams,xx,yy] = loadPotentials(sparams, 1, 'nextnano', 'xy');
xx = xx/aBstar/1E-9 + 15;
yy = yy/aBstar/1E-9;
gparams.xx = xx;
gparams.yy = yy;
gparams.nxGrid = length(xx);
gparams.nyGrid = length(yy);
sparams.nxGrid = gparams.nxGrid;
sparams.nyGrid = gparams.nyGrid;
gparams.dx = xx(2) - xx(1);
gparams.dy = yy(2) - yy(1);
[XX,YY] = meshgrid(xx,yy);
gparams.XX = XX;
gparams.YY = YY;
fprintf(1,'Done!\n');

% Find which index corresponds to where the 2DEG should be
for ii = 1:length(sparams.potentials)
    sparams.potentials(ii).pot2D = sparams.potentials(ii).pot2D/Rystar/sparams.ee/1E-3;
end

% Now we want to make the potential interpolant object (both 2D and 2DEG)
sparams.interpType = 'spline';
sparams.extrapType = 'spline';
sparams = makePotentialsInterpolantsAlt(sparams,xx,yy);

pot2DInterp = sparams.P2DInterpolant;

clearvars -except gparams pot2DInterp
%%
% Test everything loaded correctly
% voltVec = [0.3,0.3,0.2];
voltVec = [0.3,0.298925,0.2]; % 0 detuning voltage configuration for dots 1/2
% currPotential = sparams.potentials(1).pot2D;
currPotential = squeeze(pot2DInterp(getInterpolantArgument(voltVec,gparams.xx,gparams.yy)));
gparams.VV = currPotential;

aBstar = 2.16; % Convert abStar to nm
Rystar = 42.7; % Convert Ry* to meV

figure('Color','white','Position',[100,200,1800,700]);
subplot(1,3,1);
s = surf(gparams.XX,gparams.YY,squeeze(currPotential));
colormap(viridis());
set(s,'EdgeAlpha',0);
view(0,0);
xlim([min(gparams.xx),max(gparams.xx)]);
ylim([min(gparams.yy), max(gparams.yy)]);
set(gca,'TickLabelInterpreter','latex','Fontsize',14);
xlabel('x-coord [$a_B^*$]','Interpreter','latex','Fontsize',20);
zlabel('Energy [Ry*]','Interpreter','latex','Fontsize',20);

simparams;
gparams.ngridx = gparams.nxGrid;
gparams.ngridy = gparams.nyGrid;

% Get wavefunctions
[wfs, ens] = solve2DSingleElectronSE(sparams, gparams, 2);
subplot(1,3,2);
s = surf(gparams.XX,gparams.YY,wfs(:,:,1));
colormap(viridis());
set(s,'EdgeAlpha',0);
view(0,0);
xlim([min(gparams.xx),max(gparams.xx)]);
ylim([min(gparams.yy), max(gparams.yy)]);
set(gca,'TickLabelInterpreter','latex','Fontsize',14);
xlabel('x-coord [$a_B^*$]','Interpreter','latex','Fontsize',20);
zlabel('$\psi_0(x,y)$','Interpreter','latex','Fontsize',20);

subplot(1,3,3);
s = surf(gparams.XX,gparams.YY,wfs(:,:,2));
colormap(viridis());
set(s,'EdgeAlpha',0);
view(0,0);
xlim([min(gparams.xx),max(gparams.xx)]);
ylim([min(gparams.yy), max(gparams.yy)]);
set(gca,'TickLabelInterpreter','latex','Fontsize',14);
xlabel('x-coord [$a_B^*$]','Interpreter','latex','Fontsize',20);
zlabel('$\psi_1(x,y)$','Interpreter','latex','Fontsize',20);

export_fig 'ChgStab_nn++_0detuning' -m3

%%
% Find the lever arms
voltVec = [0.3,0.298925,0.2]; % 0 detuning voltage configuration for dots 1/2
% currPotential = sparams.potentials(1).pot2D;
currPotential = squeeze(pot2DInterp(getInterpolantArgument(voltVec,gparams.xx,gparams.yy)));
y0Index = getClosestArrayIndex(0,gparams.yy);
gparams.VV = currPotential;

aBstar = 2.16; % Convert abStar to nm
Rystar = 42.7; % Convert Ry* to meV

figure('Color','white','Position',[100,200,1800,700]);
subplot(1,3,1);
findpeaks(-gparams.VV(y0Index,:), gparams.xx);
set(gca,'TickLabelInterpreter','latex','Fontsize',14);
xlabel('x-coord [$a_B^*$]','Interpreter','latex','Fontsize',20);
ylabel('Energy [Ry*]','Interpreter','latex','Fontsize',20);
ylim([52,52.3]);
xlim([-40,40]);
[pks, ~] = findpeaks(-gparams.VV(y0Index,:), gparams.xx);
refDot1Pk = pks(2);
refDot2Pk = pks(3);

dV = 0.001;
voltVec = [0.3,0.298925 - dV,0.2]; % 0 detuning voltage configuration for dots 1/2
% currPotential = sparams.potentials(1).pot2D;
currPotential = squeeze(pot2DInterp(getInterpolantArgument(voltVec,gparams.xx,gparams.yy)));
gparams.VV = currPotential;

subplot(1,3,2);
findpeaks(-gparams.VV(y0Index,:), gparams.xx);
set(gca,'TickLabelInterpreter','latex','Fontsize',14);
xlabel('x-coord [$a_B^*$]','Interpreter','latex','Fontsize',20);
ylabel('Energy [Ry*]','Interpreter','latex','Fontsize',20);
ylim([52,52.3]);
xlim([-40,40]);
[pks, ~] = findpeaks(-gparams.VV(y0Index,:), gparams.xx);
dot1Pk = pks(2);
dot2Pk = pks(3);
dDot1_dV2 = -(refDot1Pk - dot1Pk)/dV;
dDot2_dV2 = -(refDot2Pk - dot2Pk)/dV;

dV = 0.001;
voltVec = [0.3 - dV,0.298925,0.2]; % 0 detuning voltage configuration for dots 1/2
% currPotential = sparams.potentials(1).pot2D;
currPotential = squeeze(pot2DInterp(getInterpolantArgument(voltVec,gparams.xx,gparams.yy)));
gparams.VV = currPotential;

subplot(1,3,3);
findpeaks(-gparams.VV(y0Index,:), gparams.xx);
set(gca,'TickLabelInterpreter','latex','Fontsize',14);
xlabel('x-coord [$a_B^*$]','Interpreter','latex','Fontsize',20);
ylabel('Energy [Ry*]','Interpreter','latex','Fontsize',20);
ylim([52,52.3]);
xlim([-40,40]);
[pks, ~] = findpeaks(-gparams.VV(y0Index,:), gparams.xx);
dot1Pk = pks(2);
dot2Pk = pks(3);
dDot1_dV1 = -(refDot1Pk - dot1Pk)/dV;
dDot2_dV1 = -(refDot2Pk - dot2Pk)/dV;

% Right now dDot_dV are in units of Ry*/V, so convert to meV/V
dDot1_dV1 = dDot1_dV1*Rystar;
dDot2_dV1 = dDot2_dV1*Rystar;
dDot1_dV2 = dDot1_dV2*Rystar;
dDot2_dV2 = dDot2_dV2*Rystar;

export_fig 'ChgStab_nn++_leverarms' -m3


%% Build double dot potential
voltVec = [0.3,0.298925,0.2]; % 0 detuning voltage configuration for dots 1/2
currPotential = squeeze(pot2DInterp(getInterpolantArgument(voltVec,gparams.xx,gparams.yy)));
y0Index = getClosestArrayIndex(0,gparams.yy);
gparams.VV = currPotential;

% sparams.dotLocations = [-dotSep/2,0;dotSep/2,0];
% [sparams.nDots,~] = size(sparams.dotLocations);

% Guess an omega value 
% omegaGuess = sqrt(Vimin/lxi^2);
omegaGuess = 0.0254;
optOmega = omegaGuess;
% optOmega = optimizeOmega(sparams,gparams,omegaGuess);
% optOmega

% Build the Harmonic orbital states
[originHOs, originOmega] = createOriginHOs(sparams,gparams,optOmega);
% plot2DBasis(gparams, originHOs, 24);

% Find the first two lowest eigenstates of the double quantum dot potential
% DQD should be at 0 detuning.
sparams.nItinerantOrbitals = 2;
[itinOrbs, itinEns] = findItinerantBasis(sparams, gparams, sparams.nItinerantOrbitals);
tc = abs(itinEns(2) - itinEns(1))/2;

% Create localized wavefunctions by taking the symmetric and antisymmetric
% combinations of two lowest eigenstates
symState = 1/sqrt(2)*(itinOrbs(1).wavefunctionMG + itinOrbs(2).wavefunctionMG);
antState = 1/sqrt(2)*(itinOrbs(1).wavefunctionMG - itinOrbs(2).wavefunctionMG);

figure('Color','white','Position',[100,200,1000,800]);
subplot(2,2,1);
surf(gparams.XX, gparams.YY, itinOrbs(1).wavefunctionMG);
set(gca,'Fontsize',14,'TickLabelInterpreter','latex');
xlabel('x-coord [$a_B^*$]','Interpreter','latex','Fontsize',20);
ylabel('y-coord [$a_B^*$]','Interpreter','latex','Fontsize',20);
xlim([min(gparams.xx),max(gparams.xx)]);
ylim([min(gparams.yy),max(gparams.yy)]);
colormap(viridis());
shading interp;
view(2);

subplot(2,2,2);
surf(gparams.XX, gparams.YY, itinOrbs(2).wavefunctionMG);
set(gca,'Fontsize',14,'TickLabelInterpreter','latex');
xlabel('x-coord [$a_B^*$]','Interpreter','latex','Fontsize',20);
ylabel('y-coord [$a_B^*$]','Interpreter','latex','Fontsize',20);
xlim([min(gparams.xx),max(gparams.xx)]);
ylim([min(gparams.yy),max(gparams.yy)]);
colormap(viridis());
shading interp;
view(2);

subplot(2,2,3);
surf(gparams.XX, gparams.YY, symState);
set(gca,'Fontsize',14,'TickLabelInterpreter','latex');
xlabel('x-coord [$a_B^*$]','Interpreter','latex','Fontsize',20);
ylabel('y-coord [$a_B^*$]','Interpreter','latex','Fontsize',20);
xlim([min(gparams.xx),max(gparams.xx)]);
ylim([min(gparams.yy),max(gparams.yy)]);
colormap(viridis());
shading interp;
view(2);

subplot(2,2,4);
surf(gparams.XX, gparams.YY, antState);
set(gca,'Fontsize',14,'TickLabelInterpreter','latex');
xlabel('x-coord [$a_B^*$]','Interpreter','latex','Fontsize',20);
ylabel('y-coord [$a_B^*$]','Interpreter','latex','Fontsize',20);
xlim([min(gparams.xx),max(gparams.xx)]);
ylim([min(gparams.yy),max(gparams.yy)]);
colormap(viridis());
shading interp;
view(2);

export_fig 'ChgStab_n++_sym_ant' -m3
%%
% Now build up a new basis of wavefunctions that we want to calculate
% coulomb interactions for
basisWFs(1).wavefunctionMG = symState;
basisWFs(2).wavefunctionMG = antState;

% Now find the transformation from HO basis to this new basis
acoeffs = findTMatrixViaInnerProd(gparams, originHOs, basisWFs);

checkBasisTransformation(sparams, gparams, originHOs, basisWFs, acoeffs);
%%
% Convert CMEs from HO basis to new basis
% load('CMEsOrigin_11_11.mat');
A = sqrt(sparams.hbar/(sparams.me*originOmega));
acoeffsTol = 1E-10;
acoeffs(acoeffs < acoeffsTol & acoeffs > -acoeffsTol) = 0;
acoeffs = sparse(acoeffs);
kronAcoeffs = kron(acoeffs,acoeffs);
CMEsNewBasis = full((kronAcoeffs*(CMEs_lib/A))*kronAcoeffs');

% Now get the CMEs we actually care about...
% <11|v|11> is dot 1's charging energy
U1 = CMEsNewBasis(1,1);
% <22|v|22> is dot 2's charging energy
U2 = CMEsNewBasis(4,4);
% <12|v|21> = <21|v|12> is the interdot charging energy
U12 = CMEsNewBasis(2,2); % Or U12 = CMEsNewBasis(3,3)

% Convert from Rydbergs to Si energy units meV
U1 = U1*Rystar; 
U2 = U2*Rystar;
U12 = U12*Rystar;

fprintf(1,'U1 = %.3E meV\n', U1);
fprintf(1,'U2 = %.3E meV\n', U2);
fprintf(1,'U12 = %.3E meV\n', U12);
fprintf(1,'tc = %.3E ueV\n', tc*Rystar/1E-3);
U1_SI = (U1*sparams.ee/1E3);
U2_SI = (U2*sparams.ee/1E3);
dot1R = 1/(U1_SI/sparams.ee^2*(8*sparams.eps0*sparams.epsR));
dot2R = 1/(U2_SI/sparams.ee^2*(8*sparams.eps0*sparams.epsR));
fprintf(1,'Dot 1 radius = %.3E nm\n',dot1R/1E-9);
fprintf(1,'Dot 2 radius = %.3E nm\n',dot2R/1E-9);

%% Simulate the charge stability

c = [0,1;0,0];
cdag = c';

c1 = kron(c,eye(2));
c1up = kron(c1,eye(4));
c1down = kron(eye(4),c1);
n1up = c1up'*c1up;
n1down = c1down'*c1down;

c2 = kron(eye(2),c);
c2up = kron(c2,eye(4));
c2down = kron(eye(4),c2);
n2up = c2up'*c2up;
n2down = c2down'*c2down;

% Charge state ordering
% (0,0) => 1
% (0,1) => 2
% (1,0) => 3
% (1,1) => 4
% (0,2) => 5
% (2,0) => 6
% (1,2) => 7
% (2,1) => 8
% (2,2) => 9

% Voltages to sweep
V1 = linspace(-0.040,0.085,100);
V2 = linspace(-0.040,0.085,100);
% V1 = -67;
% V2 = 40;
% V1 = linspace(-15.,45,100);
% V2 = linspace(-15.,45,100);
[VV1, VV2] = meshgrid(V1,V2);
colorInd = zeros(size(VV1));
virtualGates = 0;

% 
% ii = 1;
% jj = 20;

for ii = 1:length(V1)
    for jj = 1:length(V2)
            % Build current Hamiltonian
            if ~virtualGates
                mu1 = dDot1_dV1*V1(ii) + dDot1_dV2*V2(jj);
                mu2 = dDot2_dV1*V1(ii) + dDot2_dV2*V2(jj);
            else
                mu1 = V1(ii);
                mu2 = V2(jj);
            end
            
            H = -mu1*(n1up + n1down) - mu2*(n2up + n2down) +...
                U1*n1up*n1down + U2*n2up*n2down...
                + U12*(n1up*n2down + n1down*n2up + n1down*n2down + n1up*n2up);
            
            [V,D] = eig(H);
            
            % Find how many charges in each dot
            groundState = V(:,1);
%             groundState = groundState.^2; % Find the probabilty of being in each state
%             [~,ind] = max(groundState);
%             groundState = zeros(8,1);
%             groundState(ind) = 1;
            dot1Charge = groundState'*(n1up + n1down)*groundState;
            dot2Charge = groundState'*(n2up + n2down)*groundState;
            
            colorInd(ii,jj) = (dot1Charge + 1) + (dot2Charge)*3;
    end
end

figure('Color','white');
s = surf(VV1/1E-3, VV2/1E-3, colorInd);
cm = viridis(9);
colormap(cm);
set(s,'EdgeAlpha',0);
view(2);
set(gca,'Fontsize',14,'TickLabelInterpreter','latex');
xlim([min(V1),max(V1)]/1E-3);
ylim([min(V2),max(V2)]/1E-3);
if ~virtualGates
    xlabel('$V_1$ [mV]','Interpreter','latex','Fontsize',20);
    ylabel('$V_2$ [mV]','Interpreter','latex','Fontsize',20);
else
    xlabel('$V''_1$ [mV]','Interpreter','latex','Fontsize',20);
    ylabel('$V''_2$ [mV]','Interpreter','latex','Fontsize',20);
end
export_fig 'ChgStab_nn+_chgStabWithLever' -m3 

% legend{}







