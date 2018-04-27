% Always load the simulation parameters first
clear;
simparams;


% Load the potential profile to simulate
fprintf(1,'Loading potentials from file %s...\n',sparams.potFile);
[sparams,X,Y,V] = loadPotentialsFromFigure(sparams,'potentials/Brandon_1660mV.fig');


fprintf(1,'Fitting potentials to localized Harmonic Orbitals...\n');
% We want to fit the potentials to localized HOs.  We will only do this fit
% though for the lower part of the full potential (i.e. values of Vq below
% the tunnel barrier between two dots.  We only need a close fit, it does
% not have to be exact since we are just constructing a set of 1D wave
% functions in the x and y axes to use as a basis for the arbitrary single
% particle wave functions
sparams = fitDotLocationsToHarmonicWells(sparams,X,Y,V);


fprintf(1,'Finding 2D localized harmonic orbitals...\n');
sparams = solveFor2DLocalizedHOs(sparams,X,Y);

% Show the 2D LOHOs
debugHere = 0;
if sparams.verbose && debugHere
    for ii = 1:sparams.nSingleOrbitals

        if mod(ii-1,sparams.nLocalOrbitals) == 0
            fig = figure;
            pos = get(fig,'position');
            set(fig,'position',[pos(1:2)/4 pos(3)*2.5 pos(4)*2]);
        end
        currWFMG = sparams.localHOs(ii).wavefunctionMG;
        currWFNO = sparams.localHOs(ii).wavefunctionNO;
       
        subplot(2,ceil(sparams.nLocalOrbitals/2),mod(ii-1,sparams.nLocalOrbitals)+1);
        s = surf(X,Y,currWFMG);
        set(s,'edgecolor','none');
        title(sprintf('N: %d  M: %d',sparams.localHOs(ii).n,sparams.localHOs(ii).m));
        view(2);
       
        fprintf(1,'LOHO state %d norm: %0.4f \n', ii, norm(currWFNO));
        drawnow;
    end
end


fprintf(1,'Assembling the linear combination of localized harmonic orbitals...\n');
% The key item we get here is the acoeffs.  These are used in our extensive
% summation when we break down the Coulomb Matrix Element in terms of
% non-shifted HOs.  We still extract the numerically determined single
% particles as well for debugging/visualization of the code's performance.
% Now we will calculate the overlaps to change bases
% From this we get a matrix A with elements A_ij = <k_j|r_i> where k_j are
% the single particle orbitals for an arbitrary potential and r_
% are the localized HOs (i.e. shifted) and alpha_j are the HOs centered at 
% the origin.
sparams = solveLinearCombinationHOs(sparams,X,Y,V);

% Check that the acoeffs are correct by comparing to the numerical results
% A good "check" is that the energies are with 0.1 ueV
debugHere = 1;
if sparams.verbose && debugHere
    nTemp = 10;
    
    % Numerically find the 1 electron wavefunctions for the arbitrary 2D
    % potential well
    [Rwfs, Rens] = solve2DSingleElectronSE(sparams, X, Y, V, nTemp);
    
    for ii = 1:nTemp
        currWF = zeros(sparams.ngridy,sparams.ngridx);
        
        for kk = 1:sparams.nSingleOrbitals
            currWF = currWF + sparams.acoeffs(ii,kk)*sparams.localHOs(kk).wavefunctionMG;
        end

        fig = figure;
        pos = get(fig,'position');
        set(fig,'position',[pos(1:2)/4 pos(3:4)*2]);
        subplot(1,2,1);
        s = surf(X,Y,currWF);
        set(s,'edgecolor','none');
        title(sprintf('LCHO %d',ii));
        view(2);
        
        subplot(1,2,2);
        s = surf(X,Y,Rwfs(:,:,ii));
        set(s,'edgecolor','none');
        title(sprintf('Numerical %d',ii));
        view(2);
        drawnow;
        
        fprintf(1,'State %03d energies - LCHO: %0.6f eV, Numerical: %0.6f eV \n',...
            ii, sparams.linearCombinationHOs(ii).energy/sparams.ee,Rens(ii,ii)/sparams.ee);
    end
end


fprintf(1,'Calculating change of basis elements for shifted harmonic orbitals\n');
% What we are doing here is the first part of the CME calculation.  We have
% an analytical formula for the CME when we are dealing harmonic orbitals
% coming from the same well and centered at the origin.  Unfortunately, we
% have a double well system and neither of our dots are centered at the
% origin.  So we must unfold our shifted harmonic orbital wave functions in
% terms of the unshifted harmonic orbital basis.  The resulting matrix of
% bcoeffs describes this change of basis
% First, we will find our new set of basis functions centered at the
% orbital
sparams = createNonShiftedHOs(sparams,X,Y);

% Now we will calculate the overlaps to change bases
% From this we get a matrix B with elements B_ij = <r_j|alpha_i> where r_j
% are the localized HOs (i.e. shifted) and alpha_i are the HOs centered at 
% the origin.
sparams = solveShiftToNonShiftedCoeffs(sparams);

% Check that the bcoeffs are correct by comparing to the numerical results
debugHere = 1;
if sparams.verbose && debugHere
    nTemp = 10;
    
    % Numerically find the 1 electron wavefunctions for the arbitrary 2D
    % potential well
    [Rwfs, ~] = solve2DSingleElectronSE(sparams, X, Y, V, nTemp);
    
    % Build the transformation matrix to go from non-shifted HOs to single
    % particle orbitals
    transformationMatrix = sparams.acoeffs*sparams.bcoeffs;
    for jj = 1:nTemp
        currWF = zeros(sparams.ngridy,sparams.ngridx);
        for kk = 1:((sparams.maxNonShiftedHOsX+1)*(sparams.maxNonShiftedHOsY+1))
            currWF = currWF + transformationMatrix(jj,kk)*sparams.nonShiftedHOs(kk).wavefunctionMG;
        end

        fig = figure;
        pos = get(fig,'position');
        set(fig,'position',[pos(1:2)/4 pos(3:4)*2]);
        subplot(1,2,1);
        s = surf(X,Y,currWF);
        title(sprintf('A and B transformed states %d',jj));
        set(s,'edgecolor','none');
        view(2);
        
        subplot(1,2,2);
        s = surf(X,Y,Rwfs(:,:,jj));
        set(s,'edgecolor','none');
        title(sprintf('Numerical %d',jj));
        view(2);
        drawnow;

        fprintf(1,'State %03d norms - a&b coeffs: %0.4f, Numerical: %0.4f \n', jj, norm(convertMGtoNO(currWF)),...
            norm(sparams.linearCombinationHOs(jj).wavefunctionNO));
    end
end

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
%%
fprintf(1,'Constructing CMEs for single particle orbitals\n');
% Now above we only found the non shifted harmonic orbitals CMEs.  From
% those, we need to build up the actual CMEs for the single particle wave
% functions (which is what we actually need in our Hamiltonian to find J).
sparams = solveCMEsSingleParticlefromSameOrbital(sparams);
% save('exchangeTest2.mat','sparams','V','X','Y');
%%
% figure(5);
% hold on;
% % ylim([-0.1,0.04])
% % xlim([-50,50])
% ylabel('Potential [eV]','Interpreter','Latex','Fontsize',18);
% xlabel('Position [nm]','Interpreter','Latex','Fontsize',18);
% plot(X(1,:)/1E-9,potHOX2/consts.ee);
% size(wfsHOX2(:,ii))
% size(X(1,:))
% for ii = 1:(params.nLocalOrbitals)
%     plot(X(1,:)/1E-9, abs(wfsHOX2(ii,:)).^2*3/5 + engHOX2(ii,ii)/consts.ee,'Linewidth',1.5);
% end

close all;
% Plot the first 10 wavefunctions & compare to the numerical results
nTemp = 10;
[Rwfs, Rens] = solve2DSingleElectronSE(sparams, X, Y, V, nTemp);
for ii = 1:(nTemp/2)
    figure;
    tempwf = sparams.linearCombinationHOs(2*ii-1).wavefunction;
    subplot(2,2,1);
    scurr = surf(X,Y,abs(tempwf).^2);
    title(sprintf('LCHO %d',2*ii-1));
    set(scurr,'edgecolor','none');
    xlim([-50,50]*1E-9);
    ylim([-50,50]*1E-9);
    xlabel('x [m]','Interpreter','Latex','Fontsize',10);
    ylabel('y [m]','Interpreter','Latex','Fontsize',10);
    view(2);
    
    subplot(2,2,2);
    scurr = surf(X,Y,abs(Rwfs(:,:,2*ii-1).^2));
    title(sprintf('Numerical %d',2*ii-1));
    set(scurr,'edgecolor','none');
    xlim([-50,50]*1E-9);
    ylim([-50,50]*1E-9);
    xlabel('x [m]','Interpreter','Latex','Fontsize',10);
    ylabel('y [m]','Interpreter','Latex','Fontsize',10);
    view(2);

    tempwf = sparams.linearCombinationHOs(2*ii).wavefunction;
    subplot(2,2,3);
    scurr = surf(X,Y,abs(tempwf).^2);
    title(sprintf('LCHO %d',2*ii));
    set(scurr,'edgecolor','none');
    xlim([-50,50]*1E-9);
    ylim([-50,50]*1E-9);
    xlabel('x [m]','Interpreter','Latex','Fontsize',10);
    ylabel('y [m]','Interpreter','Latex','Fontsize',10);
    view(2);
    
    subplot(2,2,4);
    scurr = surf(X,Y,abs(Rwfs(:,:,2*ii).^2));
    title(sprintf('Numerical %d',2*ii));
    set(scurr,'edgecolor','none');
    xlim([-50,50]*1E-9);
    ylim([-50,50]*1E-9);
    xlabel('x [m]','Interpreter','Latex','Fontsize',10);
    ylabel('y [m]','Interpreter','Latex','Fontsize',10);
    view(2);    
    
    fprintf(1,'Est LCHO E%d = %.6d eV\n',2*ii-1,sparams.linearCombinationHOs(2*ii-1).energy/sparams.ee);
    fprintf(1,'Est Num E%d = %.6d eV\n',2*ii-1,Rens(2*ii-1,2*ii-1)/sparams.ee);  
    fprintf(1,'Est LCHO E%d = %.6d eV\n',2*ii,sparams.linearCombinationHOs(2*ii).energy/sparams.ee);
    fprintf(1,'Est Num E%d = %.6d eV\n',2*ii,Rens(2*ii,2*ii)/sparams.ee);    
end

%%
ii = 1;
figure;
title('|\Psi|^2','Interpreter','Latex');

tempwf = convertNOtoMG(wfsLCHONO(:,ii),nx,ny);
subplot(2,4,1);
scurr = surf(Xq*10^9,Yq*10^9,abs(tempwf).^2);
title(sprintf('LCHO %d',ii));
set(scurr,'edgecolor','none');
xlim([-50,50]);
ylim([-50,50]);
xlabel('x [nm]','Interpreter','Latex','Fontsize',12);
ylabel('y [nm]','Interpreter','Latex','Fontsize',12);
view(2);

subplot(2,4,2);
scurr = surf(Xq*10^9,Yq*10^9,abs(Rwfs(:,:,ii).^2));
title(sprintf('Numerical %d',ii));
set(scurr,'edgecolor','none');
xlim([-50,50]);
ylim([-50,50]);
xlabel('x [nm]','Interpreter','Latex','Fontsize',12);
ylabel('y [nm]','Interpreter','Latex','Fontsize',12);
view(2);

ii = 3;

tempwf = convertNOtoMG(wfsLCHONO(:,ii),nx,ny);
subplot(2,4,3);
scurr = surf(Xq*10^9,Yq*10^9,abs(tempwf).^2);
title(sprintf('LCHO %d',ii));
set(scurr,'edgecolor','none');
xlim([-50,50]);
ylim([-50,50]);
xlabel('x [nm]','Interpreter','Latex','Fontsize',12);
ylabel('y [nm]','Interpreter','Latex','Fontsize',12);
view(2);

subplot(2,4,4);
scurr = surf(Xq*10^9,Yq*10^9,abs(Rwfs(:,:,ii).^2));
title(sprintf('Numerical %d',ii));
set(scurr,'edgecolor','none');
xlim([-50,50]);
ylim([-50,50]);
xlabel('x [nm]','Interpreter','Latex','Fontsize',12);
ylabel('y [nm]','Interpreter','Latex','Fontsize',12);
view(2);    

ii = 9;

tempwf = convertNOtoMG(wfsLCHONO(:,ii),nx,ny);
subplot(2,4,5);
scurr = surf(Xq*10^9,Yq*10^9,abs(tempwf).^2);
title(sprintf('LCHO %d',ii));
set(scurr,'edgecolor','none');
xlim([-50,50]);
ylim([-50,50]);
xlabel('x [nm]','Interpreter','Latex','Fontsize',12);
ylabel('y [nm]','Interpreter','Latex','Fontsize',12);
view(2);

subplot(2,4,6);
scurr = surf(Xq*10^9,Yq*10^9,abs(Rwfs(:,:,ii).^2));
title(sprintf('Numerical %d',ii));
set(scurr,'edgecolor','none');
xlim([-50,50]);
ylim([-50,50]);
xlabel('x [nm]','Interpreter','Latex','Fontsize',12);
ylabel('y [nm]','Interpreter','Latex','Fontsize',12);
view(2); 

ii = 10;

tempwf = convertNOtoMG(wfsLCHONO(:,ii),nx,ny);
subplot(2,4,7);
scurr = surf(Xq*10^9,Yq*10^9,abs(tempwf).^2);
title(sprintf('LCHO %d',ii));
set(scurr,'edgecolor','none');
xlim([-50,50]);
ylim([-50,50]);
xlabel('x [nm]','Interpreter','Latex','Fontsize',12);
ylabel('y [nm]','Interpreter','Latex','Fontsize',12);
view(2); 

subplot(2,4,8);
scurr = surf(Xq*10^9,Yq*10^9,abs(Rwfs(:,:,ii).^2));
title(sprintf('Numerical %d',ii));
set(scurr,'edgecolor','none');
xlim([-50,50]);
ylim([-50,50]);
xlabel('x [nm]','Interpreter','Latex','Fontsize',12);
ylabel('y [nm]','Interpreter','Latex','Fontsize',12);
view(2); 
