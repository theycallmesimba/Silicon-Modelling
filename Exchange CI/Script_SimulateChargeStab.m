%% Load nextnano simulations

%% Build double dot potential

% Simulation parameters
simparams;

aBstar = 2.16; % Convert abStar to nm
Rystar = 42.7; % Convert Ry* to meV

% Dot potential size along x (characteristic width)
% lxi = 2.3;
lxi = 55/aBstar;
% Potential minima
Vimin = 5/Rystar;
% Dot separation (x-axis)
% dotSep = 5;
dotSep = 120/aBstar;
% Dot potential eccentricity (defined as: ly/lx)
eccentricity = 1;
% Tunnel gate potential height (wrt to Vmin)
Vtun = 14/Rystar;
% Tunnel gate width
% lxtun = 50/aBstar;
lxtun = dotSep/(2*log(2)) - lxi;

fprintf(1,'*************************\n');
fprintf(1,'**** Charge stab sim ****\n');
fprintf(1,'*************************\n');
fprintf(1,'* Dot separation = %.3f    *\n', dotSep);
fprintf(1,'* Dot width = %.3f     *\n', lxi);
fprintf(1,'* Dot ecc = %.3f       *\n', eccentricity);
fprintf(1,'* Vtun height = %.3f   *\n', Vtun);
fprintf(1,'* l_xtun = %.3f    *\n',lxtun);
fprintf(1,'*************************\n');
        
sparams.dotLocations = [-dotSep/2,0;dotSep/2,0];
[sparams.nDots,~] = size(sparams.dotLocations);
        

% Fill in grid parameters
gparams.ngridx = 300;
gparams.ngridy = 300;
% gparams.xx = linspace(-8,8,gparams.ngridx);
% gparams.yy = linspace(-8,8,gparams.ngridy);
gparams.xx = linspace(-150,150,gparams.ngridx)/aBstar;
gparams.yy = linspace(-150,150,gparams.ngridy)/aBstar;
[gparams.XX,gparams.YY] = meshgrid(gparams.xx,gparams.yy);
currVV = zeros(gparams.ngridy,gparams.ngridx);

wx = lxi;
wy = lxi*eccentricity;
dot1Pot = -Vimin*exp(-((gparams.XX - sparams.dotLocations(1,1)).^2/wx^2) -...
        (gparams.YY - sparams.dotLocations(1,2)).^2/wy^2);

dot2Pot = -Vimin*exp(-((gparams.XX - sparams.dotLocations(2,1)).^2/wx^2) -...
        (gparams.YY - sparams.dotLocations(2,2)).^2/wy^2);

wxtun = lxtun;
wytun = wy;
tunPot = Vtun*exp(-(gparams.XX).^2/wxtun^2 +...
        -(gparams.YY).^2/wytun^2);

gparams.VV = dot1Pot + dot2Pot + tunPot;

% Plot the 2D Potential at 0 bias
figure('Color','white');
surf(gparams.XX*aBstar,gparams.YY*aBstar,gparams.VV);
set(gca,'Fontsize',14,'TicklabelInterpreter','latex');
xlabel('x-axis [nm]','Fontsize',18,'Interpreter','latex');
ylabel('y-axis [nm]','Fontsize',18,'Interpreter','latex');
xlim([min(min(gparams.XX)),max(max(gparams.XX))]*aBstar);
ylim([min(min(gparams.YY)),max(max(gparams.YY))]*aBstar);
colormap(viridis());
shading interp
view(2);
title('2D Potential','Fontsize',20,'Interpreter','latex');

% Guess an omega value 
omegaGuess = sqrt(Vimin/lxi^2);
optOmega = optimizeOmega(sparams,gparams,omegaGuess);

% Build the Harmonic orbital states
[originHOs, originOmega] = createOriginHOs(sparams,gparams,optOmega);

% Find the first two lowest eigenstates of the double quantum dot potential
% DQD should be at 0 detuning.
sparams.nItinerantOrbitals = 2;
[itinOrbs, itinEns] = findItinerantBasis(sparams, gparams, sparams.nItinerantOrbitals);
tc = abs(itinEns(2) - itinEns(1))/2*Rystar;

% Create localized wavefunctions by taking the symmetric and antisymmetric
% combinations of two lowest eigenstates
symState = 1/sqrt(2)*(itinOrbs(1).wavefunctionMG + itinOrbs(2).wavefunctionMG);
antState = 1/sqrt(2)*(itinOrbs(1).wavefunctionMG - itinOrbs(2).wavefunctionMG);

figure('Color','white');
subplot(2,2,1);
surf(gparams.XX*aBstar, gparams.YY*aBstar, itinOrbs(1).wavefunctionMG);
colormap(viridis());
shading interp;
view(2);

subplot(2,2,2);
surf(gparams.XX*aBstar, gparams.YY*aBstar, itinOrbs(2).wavefunctionMG);
colormap(viridis());
shading interp;
view(2);

subplot(2,2,3);
surf(gparams.XX*aBstar, gparams.YY*aBstar, symState);
colormap(viridis());
shading interp;
view(2);

subplot(2,2,4);
surf(gparams.XX*aBstar, gparams.YY*aBstar, antState);
colormap(viridis());
shading interp;
view(2);

% Now build up a new basis of wavefunctions that we want to calculate
% coulomb interactions for
basisWFs(1).wavefunctionMG = symState;
basisWFs(2).wavefunctionMG = antState;

% Now find the transformation from HO basis to this new basis
acoeffs = findTMatrixViaInnerProd(gparams, originHOs, basisWFs);

checkBasisTransformation(sparams, gparams, originHOs, basisWFs, acoeffs);

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
fprintf(1,'tc = %.3E meV\n', tc);

% plot(gparams.XX(floor(end/2),:),symState(floor(end/2),:));
[pks,locs,w,p] = findpeaks(symState(floor(end/2),:), gparams.XX(floor(end/2),:));
fprintf(1,'Dot width = %.3f\n',w);
%% Simulate the charge stability

c = [0,1;0,0];
cdag = c';

c1 = kron(c,eye(2));
c1up = kron(c1,eye(2));
c1down = kron(eye(2),c1);
n1up = c1up'*c1up;
n1down = c1down'*c1down;

c2 = kron(eye(2),c);
c2up = kron(c2,eye(2));
c2down = kron(eye(2),c2);
n2up = c2up'*c2up;
n2down = c2down'*c2down;

cm = viridis(9);
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
V1 = linspace(-5,25,100);
V2 = linspace(-5,25,100);
[VV1, VV2] = meshgrid(V1,V2);
colorInd = zeros(size(VV1));
% 
% ii = 1;
% jj = 20;

for ii = 1:length(V1)
    for jj = 1:length(V2)
            % Build current Hamiltonian
            H = -V1(ii)*(n1up + n1down) - V2(jj)*(n2up + n2down) +...
                U1*n1up*n1down + U2*n2up*n2down + U12*(n1up*n2down +...
                n1up*n2down + n1down*n2down + n1down*n2up);
            
            [V,D] = eig(H);
            
            % Find how many charges in each dot
            groundState = V(:,1);
            dot1Charge = groundState'*(n1up + n1down)*groundState;
            dot2Charge = groundState'*(n2up + n2down)*groundState;
            
            colorInd(ii,jj) = (dot1Charge + 1) + (dot2Charge)*2;
    end
end

s = surf(VV1, VV2, colorInd);
cm = viridis(9);
colormap(cm);
set(s,'EdgeAlpha',0);
view(2);







