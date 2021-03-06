Initialize_Potential;

sparams.maxOriginHOsX = 14;
sparams.maxOriginHOsY = 14;

% Simulation parameters
% Dot potential size along x (characteristic width)
lxi = 2.3;
% Potential minima
Vimin = 20;
% Dot separation (x-axis)
dotSep = 4;
% Number of canonical orbitals to include
nItinOrbs = 10;
% Dot potential eccentricity (defined as: ly/lx)
eccentricity = 1.5;
% Dot bias (bias applied to right dot)
bias = 0;
% Tunnel gate potential height (wrt to left dot potential minimum)
Vtun = 0;
% Tunnel gate width
lxtunFrac = 1/2;

fprintf(1,'*************************\n');
fprintf(1,'* Dot separation = %d    *\n', dotSep);
fprintf(1,'* Dot width = %.3f     *\n', lxi);
fprintf(1,'* N canonical orbs = %d *\n', nItinOrbs);
fprintf(1,'* Dot ecc = %.3f       *\n', eccentricity);
fprintf(1,'* Dot bias = %.3f     *\n', bias);
fprintf(1,'* Vtun height = %.3f   *\n', Vtun);
fprintf(1,'* l_xtun/l_x = %.3f    *\n',lxtunFrac);
fprintf(1,'*************************\n');

sparams.dotLocations = [-dotSep/2,0;dotSep/2,0];
% sparams.dotLocations = [-dotSep/2,0;dotSep/2,0;-dotSep*3/2,0;dotSep*3/2,0];
[sparams.nDots,~] = size(sparams.dotLocations);

omega = sqrt(Vimin/lxi^2);

% Fill in grid parameters
gparams.ngridx = 200;
gparams.ngridy = 200;
gparams.xx = linspace(-8,8,gparams.ngridx)*effaB;
gparams.yy = linspace(-8,8,gparams.ngridy)*effaB;
[gparams.XX,gparams.YY] = meshgrid(gparams.xx,gparams.yy);
currVV = zeros(gparams.ngridy,gparams.ngridx);

wx = lxi;
wy = lxi*eccentricity;
dot1Pot = (-Vimin-bias/2)*exp(-((gparams.XX - sparams.dotLocations(1,1)).^2/wx^2) -...
    (gparams.YY - sparams.dotLocations(1,2)).^2/wy^2);

dot2Pot = (-Vimin+bias/2)*exp(-((gparams.XX - sparams.dotLocations(2,1)).^2/wx^2) -...
    (gparams.YY - sparams.dotLocations(2,2)).^2/wy^2);

wxtun = wx*lxtunFrac;
wytun = wy;
tunPot = Vtun*exp(-(gparams.XX).^2/wxtun^2 +...
    -(gparams.YY).^2/wytun^2);

if sparams.nDots > 2
    dot1Pot = dot1Pot + (-Vimin)*exp(-((gparams.XX - sparams.dotLocations(3,1)).^2/wx^2) -...
        (gparams.YY - sparams.dotLocations(3,2)).^2/wy^2);
    
    tunPot = tunPot + Vtun*exp(-(gparams.XX - dotSep).^2/wxtun^2 +...
    -(gparams.YY).^2/wytun^2);

    tunPot = tunPot + Vtun*exp(-(gparams.XX + dotSep).^2/wxtun^2 +...
    -(gparams.YY).^2/wytun^2);

    dot1Pot = dot1Pot + (-Vimin)*exp(-((gparams.XX - sparams.dotLocations(4,1)).^2/wx^2) -...
        (gparams.YY - sparams.dotLocations(4,2)).^2/wy^2);
end

gparams.VV = dot1Pot + dot2Pot + tunPot;

plotMeshgrid(gparams, gparams.VV);

sparams.nItinerantOrbitals = nItinOrbs;

sparams.nOriginHOs = sparams.maxOriginHOsX*sparams.maxOriginHOsY;

omegaGuess = abs(mean(sparams.fittedPotentialParameters(:,1)));
fprintf(1,'Optimizing origin harmonic orbital omega...  ');
optOmega = optimizeOmega(sparams,gparams,omegaGuess);
fprintf(1,'Done!\n\n');

%**************************************%
fprintf(1,'Finding 2D harmonic orbitals at origin...  ');
% Create a new basis of HOs centered at the origin of our dots
[originHOs, originOmega] = createOriginHOs(sparams,gparams,optOmega);
fprintf(1,'Done!\n\n');

%**************************************%
fprintf(1,'Finding A matrix...  ');
nStates = sparams.nOriginHOs;
basisToUse = originHOs;

[itinOrbs, itinEns] = findItinerantBasis(sparams, gparams, sparams.nItinerantOrbitals);
acoeffs = findTMatrixViaInnerProd(gparams, basisToUse, itinOrbs);

% Now find the itinerant basis energies using acoeffs
sparams.LCHOEnergies = zeros(1,sparams.nItinerantOrbitals);
full2DLap = make2DSELap(sparams,gparams);
overlap = zeros(1,sparams.nItinerantOrbitals);
for ii = 1:sparams.nItinerantOrbitals
    tempwf = zeros(gparams.ngridy*gparams.ngridx,1);

    for jj = 1:nStates
        tempwf = tempwf + acoeffs(ii,jj)*originHOs(jj).wavefunctionNO;
    end
    
    overlap(ii) = getInnerProduct2D(itinOrbs(ii).wavefunctionMG, ...
        convertNOtoMG(tempwf,gparams.ngridx,gparams.ngridy), gparams.XX, gparams.YY);

    sparams.LCHOEnergies(ii) = getInnerProduct2D(itinOrbs(ii).wavefunctionMG,...
        convertNOtoMG(full2DLap*tempwf,gparams.ngridx,gparams.ngridy), gparams.XX, gparams.YY);
end

% Truncate A in accordance with how many itinerant orbitals we want
acoeffs = acoeffs(1:sparams.nItinerantOrbitals,:);
sparams.LCHOEnergies = sparams.LCHOEnergies(1:sparams.nItinerantOrbitals);
LCHOEns = sparams.LCHOEnergies;

% Check the change of basis between the Horigin states and the itinerant
% basis
nStatesToCompare = sparams.nItinerantOrbitals;
basisToT = originHOs;
basisToCompare = findItinerantBasis(sparams, gparams, nStatesToCompare);

checkBasisTransformation(sparams, gparams, basisToT, basisToCompare, acoeffs);
fprintf(1,'Done!\n\n');

fprintf(1,'Closeness of itinerant orbital states.\n');
for ii = 1:sparams.nItinerantOrbitals
    fprintf(1,'State %d: 1-<\x03be''|\x03be> = %.3E\n',ii,1-overlap(ii));
    fprintf(1,'\x03b5 = %.3E, \x03b5'' = %.3E, %%E diff = %.3E\n',...
        itinEns(ii), sparams.LCHOEnergies(ii), (itinEns(ii) - sparams.LCHOEnergies(ii))/itinEns(ii));
end





