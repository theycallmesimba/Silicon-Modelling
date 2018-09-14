RytoJ = 5.93E-3*1.602E-19; % [1 Ry = 5.93 meV]
aBtom = 97.9E-10; % [1 aB* = 97.9 A]
RytoJ = 1;
aBtom = 1;
simparams;
Vimin = 10;
Vimin = Vimin*RytoJ;
di = 2.3*aBtom;
if RytoJ == 1 && aBtom == 1
    omega = sqrt(Vimin/di^2);
else
    omega = sqrt(2*Vimin/(sparams.me*di^2));
end

sparams.ngridx = 150;
sparams.ngridy = 150;
xx = linspace(-7.5,7.5,sparams.ngridx)*aBtom;
yy = linspace(-7.5,7.5,sparams.ngridy)*aBtom;
[XX,YY] = meshgrid(xx,yy);

VV = zeros(sparams.ngridx,sparams.ngridy);
for ii = 1:sparams.nDots
    VV = VV + -Vimin*exp(-((XX - sparams.dotLocations(ii,1)).^2 + (YY - sparams.dotLocations(ii,2)).^2)/di^2);
    sparams.fittedPotentialParameters(ii,:) = [omega, Vimin, sparams.dotLocations(ii,:)];
end

% figure;
% s = surf(XX,YY,VV);
% set(s,'Edgecolor','none');

fprintf(1,'Finding 2D localized harmonic orbitals...\n');
sparams = solveFor2DLocalizedHOs(sparams,XX,YY);

% Show the 2D LOHOs and check that they are properly normalized
debugHere = 1;
if sparams.verbose && debugHere
    checkLOHOStates(sparams,XX,YY);
end
%%
fprintf(1,'Performing Loewdin Orthonormalization...\n');
% We want to obtain a transformation matrix that rotates from a basis
% composed of localized harmonic orbitals and our iterant basis of orbitals
% for the arbitrary potential.  If we just take the localized HO basis of a
% single dot, it is indeed orthonormal.  But once we have more than one
% dot, composing a basis of many of these localized HO basis sets will not
% be orthonormal as the wave functions will have overlap.  To get around
% this, we do a Loewdin orthonormalization technique.  The first part of
% this is to get an overlap matrix S from which we will transform our basis
% set {\ket{a_new}} = S^{1/2}{\ket{a}}
sparams = solveLoewdinOrthonormalization(sparams, XX, YY);

% Check orthonormality of the states
debugHere = 1;
if sparams.verbose && debugHere
    checkLoeLOHOStates(sparams,XX,YY)
end
%%
[sparams, HLOHO] = solveLinearCombinationHOs(sparams,XX,YY,VV);

% Check what the LCHO orbitals look like and give their energies
debugHere = 1;
if sparams.verbose && debugHere
    checkLCHOStates(sparams, XX, YY);
end
%%

% Create a new basis of HOs centered at the origin of our dots
sparams = createOriginHOs(sparams,XX,YY);

% Show the 2D nonShiftLOHOs and check that they are properly normalized
debugHere = 1;
if sparams.verbose && debugHere
    checkOriginHOStates(sparams,XX,YY);
end
%%

sparams = solveLoeToOriginCoeffs(sparams,XX,YY);

% Check that the bcoeffs are correct by trying to build the Loe states
debugHere = 1;
if sparams.verbose && debugHere
    checkBMatrix(sparams,XX,YY);
end
%%
fprintf(1,'Checking conversion between basis sets...\n');
full2DLap = make2DSELap(sparams,XX,YY,VV);

Horigin = zeros(sparams.nOriginHOs);
for ii = 1:sparams.nOriginHOs
    for jj = 1:sparams.nOriginHOs    
        currWFLeft = sparams.originHOs(ii).wavefunctionMG;
        currWFRight = sparams.originHOs(jj).wavefunctionNO;
        
        Horigin(ii,jj) = getInnerProduct2D(currWFLeft,...
                convertNOtoMG(full2DLap*currWFRight,sparams.ngridx,sparams.ngridy),XX,YY);
    end
end

% HLOHO = zeros(sparams.nSingleOrbitals);
% for ii = 1:sparams.nSingleOrbitals
%     for jj = 1:sparams.nSingleOrbitals
%         currWFLeft = sparams.localHOs(ii).wavefunctionMG;
%         currWFRight = sparams.localHOs(jj).wavefunctionNO;
%         
%         HLOHO(ii,jj) = getInnerProduct(currWFLeft,...
%                 convertNOtoMG(full2DLap*currWFRight,sparams.ngridx,sparams.ngridy),XX,YY);
%     end
% end

[~,Hactual] = solve2DSingleElectronSE( sparams, XX, YY, VV, 10);

HLOHOA = sparams.acoeffs*HLOHO*sparams.acoeffs';
[~,perm] = sort(diag(HLOHOA));
HLOHOA = HLOHOA(perm,perm);

HoriginB = sparams.bcoeffs*Horigin*sparams.bcoeffs';
HoriginAB = sparams.acoeffs*HoriginB*sparams.acoeffs';
[~,perm] = sort(diag(HoriginAB));
HoriginAB = HoriginAB(perm,perm);

fprintf(1,'Done!\n');
%%
fprintf(1,'Evaluating CMEs for non shifted harmonic orbitals\n');
% Now we will take all of the non shifted harmonic orbitals and evaulte the
% CMEs for all of their interactions.  A lot of these should end up being
% 0, but this will be the most computationally intensive part of the
% calculation.
tic;
sparams = solveCMEsSameOrbital(sparams);
% Save because that took a very long time
% save('exchangeTest1.mat','sparams','V','X','Y');
toc;



