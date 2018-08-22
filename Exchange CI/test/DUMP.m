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
debugHere = 0;
if sparams.verbose && debugHere
    for ii = 1:sparams.nSingleOrbitals
        currWFMG = sparams.localHOs(ii).wavefunctionMG;
       
        if mod(ii-1,8) == 0
            fig = figure;
            pos = get(fig,'position');
            set(fig,'position',[pos(1:2)/4 pos(3)*1.8 pos(4)*1.3]);
        end
        subplot(2,4,mod(ii-1,8)+1);
        
        s = surf(XX,YY,currWFMG);
        set(s,'edgecolor','none');
        title(sprintf('Dot: %d  N: %d  M: %d',sparams.localHOs(ii).dot,...
            sparams.localHOs(ii).n,sparams.localHOs(ii).m));
        view(2);
        colormap(jet);
       
        fprintf(1,'LOHO state %03d - dot: %d, n: %d, m: %d, energy: %0.3E [J], norm: %0.3E \n',...
            ii, sparams.localHOs(ii).dot, sparams.localHOs(ii).n,sparams.localHOs(ii).m,...
            sparams.localHOs(ii).energy, getInnerProduct(currWFMG,currWFMG,XX,YY));
        drawnow;
    end
end

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
debugHere = 0;
if sparams.verbose && debugHere
    for ii = 1:sparams.nSingleOrbitals
        for jj = 1:sparams.nSingleOrbitals
            printNorm = 0;
            currNorm = getInnerProduct(sparams.sLocalHOs(ii).wavefunctionMG,...
                sparams.sLocalHOs(jj).wavefunctionMG,X,Y);
            % We will only print out the norm if it is below our acceptable
            % norm tolerance threshold.  Since we are doing things
            % numerically there is bound to be some error in our
            % orthogonality not being explicitly 0
            if ii == jj && abs(currNorm - 1) >= sparams.normThreshold
                printNorm = 1;
                fprintf(1,'Norm of Loe states not below threshold <%d|%d> = %g\n',ii,jj,currNorm);
            elseif ii ~= jj && abs(currNorm - 0) >= sparams.normThreshold
                fprintf(1,'Norm of Loe states not below threshold <%d|%d> = %g\n',ii,jj,currNorm);
            end   
        end
        if mod(ii-1,10) == 0
            fig = figure;
            pos = get(fig,'position');
            set(fig,'position',[pos(1:2)/4 pos(3)*2.5 pos(4)*2]);
        end
       
        subplot(2,5,mod(ii-1,10)+1);
        s = surf(X,Y,sparams.sLocalHOs(ii).wavefunctionMG);
        set(s,'edgecolor','none');
        title(sprintf('%Leo LOHOs: %d',ii));
        view(2);
        colormap(hot);
        drawnow;
    end
end

[sparams, HLOHO, SLOHO] = solveLinearCombinationHOs(sparams,XX,YY,VV);

% Check what the LCHO orbitals look like and give their energies
debugHere = 0;
if sparams.verbose && debugHere

    for ii = 1:sparams.nSingleOrbitals     
        currWF = sparams.linearCombinationSHOs(ii).wavefunctionMG;

        if mod(ii-1,8) == 0
            fig = figure;
            pos = get(fig,'position');
            set(fig,'position',[pos(1:2)/4 pos(3)*1.8 pos(4)*1.3]);
        end
        subplot(2,4,mod(ii-1,8)+1);
        s = surf(XX,YY,currWF);
        set(s,'edgecolor','none');
        title(sprintf('LCHO %d',ii));
        view(2);
        colormap(jet);

        fprintf(1,'LCHO state %03d - energy: %0.3E [J]\n',...
            ii, sparams.linearCombinationSHOs(ii).energy);
    end
end


sparams = createNonShiftedHOs(sparams,XX,YY);

% Show the 2D nonShiftLOHOs and check that they are properly normalized
debugHere = 0;
if sparams.verbose && debugHere
    for ii = 1:sparams.nNonShiftedHOs
        currWFMG = sparams.nonShiftedHOs(ii).wavefunctionMG;
       
        if mod(ii-1,8) == 0
            fig = figure;
            pos = get(fig,'position');
            set(fig,'position',[pos(1:2)/4 pos(3)*1.8 pos(4)*1.3]);
        end
        subplot(2,4,mod(ii-1,8)+1);
        
        s = surf(XX,YY,currWFMG);
        set(s,'edgecolor','none');
        title(sprintf('Non-shifted HOs N: %d  M: %d',...
            sparams.nonShiftedHOs(ii).n,sparams.nonShiftedHOs(ii).m));
        view(2);
        colormap(jet);
       
        fprintf(1,'Non-shifted HO state %03d - n: %d, m: %d, energy: %0.3E [J], norm: %0.3E \n',...
            ii, sparams.nonShiftedHOs(ii).n,sparams.nonShiftedHOs(ii).m,...
            sparams.nonShiftedHOs(ii).energy, getInnerProduct(currWFMG,currWFMG,XX,YY));
        drawnow;
    end
end


sparams = solveShiftToNonShiftedCoeffs(sparams,XX,YY);

% Check that the bcoeffs are correct by trying to build the Loe states
debugHere = 0;
if sparams.verbose && debugHere
    for jj = 1:sparams.nSingleOrbitals
        currWF = zeros(sparams.ngridy,sparams.ngridx);
        for kk = 1:sparams.nNonShiftedHOs
            currWF = currWF + sparams.bcoeffs(jj,kk)*sparams.nonShiftedHOs(kk).wavefunctionMG;
        end

        fig = figure;
        pos = get(fig,'position');
        set(fig,'position',[pos(1:2)/4 pos(3)*1.8 pos(4)*1.3]);
        subplot(1,2,1);
        s = surf(XX,YY,currWF);
        title(sprintf('B transformed %d',jj));
        set(s,'edgecolor','none');
        colormap(hot)
        view(2);
        
        subplot(1,2,2);
        s = surf(XX,YY,sparams.sLocalHOs(jj).wavefunctionMG);
        set(s,'edgecolor','none');
        title(sprintf('LOHO %d',jj));
        view(2);
        colormap(hot);
        drawnow;

        fprintf(1,'State %03d norms - b transformed: %0.4f, Loe LOHO states: %0.4f \n',...
            jj, getInnerProduct(currWF,currWF,XX,YY),...
            getInnerProduct(sparams.linearCombinationSHOs(jj).wavefunctionMG,...
            sparams.linearCombinationSHOs(jj).wavefunctionMG,XX,YY));
    end
end

fprintf(1,'Checking conversion between basis sets...\n');
full2DLap = make2DSELap(sparams,XX,YY,VV);

HnonShift = zeros(sparams.nNonShiftedHOs);
for ii = 1:sparams.nNonShiftedHOs
    for jj = 1:sparams.nNonShiftedHOs    
        currWFLeft = sparams.nonShiftedHOs(ii).wavefunctionMG;
        currWFRight = sparams.nonShiftedHOs(jj).wavefunctionNO;
        
        HnonShift(ii,jj) = getInnerProduct(currWFLeft,...
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

HnonShiftedB = sparams.bcoeffs*HnonShift*sparams.bcoeffs';
HnonShiftedAB = sparams.acoeffs*HnonShiftedB*sparams.acoeffs';
[~,perm] = sort(diag(HnonShiftedAB));
HnonShiftedAB = HnonShiftedAB(perm,perm);


fprintf(1,'DONE!\n');



