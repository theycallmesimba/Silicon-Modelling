% Always load the simulation parameters first
% clear;
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

% Show the 2D LOHOs and check that they are properly normalized
debugHere = 1;
if sparams.verbose && debugHere
    for ii = 1:sparams.nSingleOrbitals

        if mod(ii-1,sparams.nLocalOrbitals) == 0
            fig = figure;
            pos = get(fig,'position');
            set(fig,'position',[pos(1:2)/4 pos(3)*2.5 pos(4)*2]);
        end
        currWFMG = sparams.localHOs(ii).wavefunctionMG;
       
        subplot(2,ceil(sparams.nLocalOrbitals/2),mod(ii-1,sparams.nLocalOrbitals)+1);
        s = surf(X,Y,currWFMG);
        set(s,'edgecolor','none');
        title(sprintf('N: %d  M: %d',sparams.localHOs(ii).n,sparams.localHOs(ii).m));
        view(2);
        colormap(jet);
       
        fprintf(1,'LOHO state %d norm: %0.4f \n', ii,getInnerProduct(currWFMG,currWFMG,X,Y));
        drawnow;
    end
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
sparams = solveLoewdinOrthonormalization(sparams, X, Y);

% Check orthonormality of the states
debugHere = 1;
if sparams.verbose && debugHere
    for ii = 1:sparams.nSingleOrbitals
        for jj = 1:sparams.nSingleOrbitals
            printNorm = 0;
            currNorm = getInnerProduct(sparams.sLocalHOs(ii).wavefunctionMG,...
                sparams.sLocalHOs(jj).wavefunctionMG,X,Y);
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
%%

fprintf(1,'Assembling the linear combination of localized harmonic orbitals...\n');
% The key item we get here is the acoeffs.  These are used in our extensive
% summation when we break down the Coulomb Matrix Element in terms of
% non-shifted HOs.  We still extract the numerically determined single
% particles as well for debugging/visualization of the code's performance.
% Now we will calculate the overlaps to change bases
% From this we get a matrix A with elements A_ij = <r_j|k_i> where k are
% the single particle orbitals for an arbitrary potential and r
% are the Loewdin transformed localized HOs
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
        currWF = sparams.linearCombinationSHOs(ii).wavefunctionMG;
%         currWF = sparams.linearCombinationNSHOs(ii).wavefunctionMG;

        fig = figure;
        pos = get(fig,'position');
        set(fig,'position',[pos(1:2)/4 pos(3:4)*2]);
        subplot(1,2,1);
        s = surf(X,Y,currWF);
        set(s,'edgecolor','none');
        title(sprintf('LCHO %d',ii));
        view(2);
        colormap(jet);
        
        subplot(1,2,2);
        s = surf(X,Y,Rwfs(:,:,ii));
        set(s,'edgecolor','none');
        title(sprintf('Numerical %d',ii));
        view(2);
        colormap(jet);
        drawnow;
        
        fprintf(1,'State %03d energies - LCHO: %0.6f eV, Numerical: %0.6f eV \n',...
            ii, sparams.linearCombinationSHOs(ii).energy/sparams.ee,Rens(ii,ii)/sparams.ee);
%         fprintf(1,'State %03d energies - LCHO: %0.5f eV, Numerical: %0.5f eV \n',...
%             ii, sparams.linearCombinationNSHOs(ii).energy/sparams.ee,Rens(ii,ii)/sparams.ee);
    end
end
%%

fprintf(1,'Building basis of non shifted harmonic orbitals...\n');
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

% Check orthonormality of the states
debugHere = 1;
if sparams.verbose && debugHere
    % Will only print out orthogonality conditions if they are not below
    % threshold
    for ii = 1:sparams.nNonShiftedHOs
        for jj = 1:sparams.nNonShiftedHOs
            currNorm = getInnerProduct(sparams.nonShiftedHOs(ii).wavefunctionMG,...
                sparams.nonShiftedHOs(jj).wavefunctionMG,X,Y);
            if ii == jj && abs(currNorm - 1) >= sparams.normThreshold
                fprintf(1,'Norm of non-shifted HO states not below threshold <%d|%d> = %g\n',ii,jj,currNorm);
            elseif ii ~= jj && abs(currNorm - 0) >= sparams.normThreshold
                fprintf(1,'Norm of non-shifted HO states not below threshold <%d|%d> = %g\n',ii,jj,currNorm);
            end   
        end
        if mod(ii-1,10) == 0
            fig = figure;
            pos = get(fig,'position');
            set(fig,'position',[pos(1:2)/4 pos(3)*2.5 pos(4)*2]);
        end
       
        subplot(2,5,mod(ii-1,10)+1);
        s = surf(X,Y,sparams.nonShiftedHOs(ii).wavefunctionMG);
        set(s,'edgecolor','none');
        title(sprintf('%Non-Shifted LOHOs: %d',ii));
        view(2);
        colormap(hot);
        drawnow;
    end
end
%%

fprintf('Calculating change of basis elements for shifted harmonic orbitals...\n')
% Now we will calculate the overlaps to change bases
% From this we get a matrix B with elements B_ij = <alpha_j|r_i> where r
% are the localized HOs (i.e. shifted) and alpha are the HOs centered at 
% the origin.
sparams = solveShiftToNonShiftedCoeffs(sparams,X,Y);

% Check that the bcoeffs are correct by trying to build the Loe states
debugHere = 1;
if sparams.verbose && debugHere
    nTemp = 12;
    
    for jj = 1:nTemp
        currWF = zeros(sparams.ngridy,sparams.ngridx);
        for kk = 1:sparams.nNonShiftedHOs
            currWF = currWF + sparams.bcoeffs(jj,kk)*sparams.nonShiftedHOs(kk).wavefunctionMG;
        end

        fig = figure;
        pos = get(fig,'position');
        set(fig,'position',[pos(1:2)/4 pos(3:4)*2]);
        subplot(1,2,1);
        s = surf(X,Y,currWF);
        title(sprintf('B transformed %d',jj));
        set(s,'edgecolor','none');
        colormap(hot)
        view(2);
        
        subplot(1,2,2);
        s = surf(X,Y,sparams.sLocalHOs(jj).wavefunctionMG);
        set(s,'edgecolor','none');
        title(sprintf('Loe LOHO %d',jj));
        view(2);
        colormap(hot);
        drawnow;

        fprintf(1,'State %03d norms - b transformed: %0.4f, Loe LOHO states: %0.4f \n',...
            jj, getInnerProduct(currWF,currWF,X,Y),...
            getInnerProduct(sparams.linearCombinationSHOs(jj).wavefunctionMG,...
            sparams.linearCombinationSHOs(jj).wavefunctionMG,X,Y));
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
