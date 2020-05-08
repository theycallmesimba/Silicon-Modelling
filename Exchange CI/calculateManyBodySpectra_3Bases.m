function [manyBody_evecs, manyBody_ens] = calculateManyBodySpectra_3Bases(...
    sparams, gparams, optOmegaFlag, debug)
%CALCULATEMANYBODYSPECTRA This function calculates the many electron energy spectra given
%an arbitrary potential landscape.  The energy spectra is found using a
%LCHO-CI approach where we approximate the itinerant wavefunctions of the
%potential using a basis of Harmonic orbitals centered at the origin.  This
%calculation will load a pre-calculated set of Coulomb Matrix Elements if
%it is available (if not, it will calculate them but this computation time
%is costly).
%   Input parameters:
%   sparams: structure containing all of the simulation parameters.  Is
%   setup using the simparams.m file.
%   gparams: structure containing all relevant information for the grid.
%   TODO: add more detail on gparams variables
%   debug: flag either 0 or 1 used to initiate debug mode or not

% NOTE: RIGHT NOW THE CODE IS USING 3 BASES TO CALCULATE ENERGIES TO
% COMPARE WITH MAREK'S PAPER.  THIS SHOULD CHANCE TO JUST 2 BASES AFTER
% TESTING IS COMPLETE.
    fprintf(1,'Begining many body energy calculation\n\n');

    %**************************************%
    fprintf(1,'Finding 2D localized harmonic orbitals...  ');
    LOHOs = solveFor2DLocalizedHOs(sparams,gparams);

    % Show the 2D LOHOs and check that they are properly normalized (should NOT
    % be orthogonal)
    if sparams.verbose && debug
        basisToCheck = LOHOs;
        nStates = sparams.nSingleOrbitals;
        plot2DBasis(gparams,basisToCheck,nStates);
        check2DBasisOrthogonality(sparams,gparams,basisToCheck,nStates);
        check2DBasisNormality(sparams,gparams,basisToCheck,nStates);
    end
    fprintf(1,'Done!\n\n');
    
    %**************************************%
    fprintf(1,'Finding A matrix...  ');
    nStates = sparams.nSingleOrbitals;
    basisToUse = LOHOs;
    [acoeffs, ~, sparams.LCHOEnergies] = findTMatrixViaHamiltonian(sparams, gparams, basisToUse, nStates);

    % Check if the transformation is orthonormal
    if sparams.verbose && debug
        checkTMatrix(sparams, acoeffs);
    end
    % Check the change of basis between the HLOHO states and the itinerant
    % basis
    if sparams.verbose && debug
        nStatesToCompare = sparams.nSingleOrbitals;
        basisToT = LOHOs;
        basisToCompare = findItinerantBasis(sparams, gparams, nStatesToCompare);

        checkBasisTransformation(sparams, gparams, basisToT, basisToCompare, acoeffs);
    end
    fprintf(1,'Done!\n\n');

    %**************************************%
    omegaGuess = abs(mean(sparams.fittedPotentialParameters(:,1)));
    if optOmegaFlag
        fprintf(1,'Optimizing origin harmonic orbital omega...  ');
        optOmega = optimizeOmega(sparams,gparams,omegaGuess);
        fprintf(1,'Done!\n\n');
    else
        optOmega = omegaGuess;
    end
    
    %**************************************%
    fprintf(1,'Finding 2D harmonic orbitals at origin...  ');
    % Create a new basis of HOs centered at the origin of our dots
    [originHOs, originOmega] = createOriginHOs(sparams,gparams,optOmega);

    % Check orthonormality of the states
    if sparams.verbose && debug
        basisToCheck = originHOs;
        nStates = sparams.nOriginHOs;
        plot2DBasis(gparams,basisToCheck,16);
        check2DBasisOrthogonality(sparams,gparams,basisToCheck,nStates);
        check2DBasisNormality(sparams,gparams,basisToCheck,nStates);
    end
    fprintf(1,'Done!\n\n');

    %**************************************%
    fprintf(1,'Finding B matrix...  ');
    basisInit = originHOs;
    basisEnd = LOHOs;
    bcoeffs = findTMatrixViaInnerProd(gparams, basisInit, basisEnd);

    % Check if the transformation is orthonormal
    debug = 0;
    if sparams.verbose && debug
        checkTMatrix(sparams, bcoeffs);
    end
    % Check the change of basis between the Horigin states and the HloeLOHO
    % basis
    if sparams.verbose && debug
        basisToT = originHOs;
        basisToCompare = LOHOs;
        checkBasisTransformation(sparams, gparams, basisToT, basisToCompare, bcoeffs);
    end
    % Check the change of basis between the Horigin states and the itinerant
    % basis
    if sparams.verbose && debug
        nStatesToCompare = sparams.nSingleOrbitals;
        basisToT = originHOs;
        basisToCompare = findItinerantBasis(sparams, gparams, nStatesToCompare);
        checkBasisTransformation(sparams, gparams, basisToT, basisToCompare, acoeffs*bcoeffs);
    end
    fprintf(1,'Done!\n\n');

    %**************************************%
    if exist(sparams.CMEsLib_fPath, 'file') == 2
        % CME library exists.
        fprintf(1,'Found the CMEs_lib file!\n');
         if exist('CMEs_lib', 'var') == 1
            % Library already loaded
         else
            % Library needs to be loaded
            fprintf(1,'Loading CME library for calculation...  ');
            load(sparams.CMEsLib_fPath);
            assignin('base','CMEs_lib',CMEs_lib);
         end
    else
        % CME library does not exist.
        fprintf(1,'Could not find CMEs_lib file!\n');
        fprintf(1,'Evaluating CMEs for non shifted harmonic orbitals...\n');
        % Now we will take all of the non shifted harmonic orbitals and 
        % evaluate the CMEs for all of their interactions.  A lot of these 
        % should end up being 0 (75% actually), but this is the most 
        % computationally intensive part of the calculation.
        tic;
        CMEs_lib = constructOriginCMEsLib(sparams);
        elapsedTime = toc;
        fprintf(1,'Evaluating the origin CMEs took %.4f sec\n', elapsedTime);
        
        % Save because that took a very long time
        saveString = sprintf('CMEsOrigin_%d_%d_%s',CMEs_lib, ...
            sparams.maxOriginHOsX, sparams.maxOriginHOsY,...
            datestr(datetime('now'),'yy-mm-dd'));
        nStatesX_CMEsOrigin = sparams.maxOriginHOsX;
        nStatesY_CMEsOrigin = sparams.maxOriginHOsY;
        save(saveString,'CMEs_lib',nStatesX_CMEsOrigin,nStatesY_CMEsOrigin);
    end
    fprintf(1,'Done!\n\n');
    
    % Get a subset of the CMEs if the given nOrigins parameter in the
    % simparams file is less than what we've solved for in the library
    % already. Useful for checking convergence wrt number of orbitals. Note
    % that since it was defined in simparams, that the B matrix will
    % already be the correct dimensionality for the subset of CMEs_lib.
    CMEs_lib_sub = getSubsetOfCMEs(sparams, CMEs_lib,...
        nStatesX_CMEsOrigin*nStatesY_CMEsOrigin);
    
    % Scale the CMEs to match the origin HOs used to form transformation
    % matrix B
    A = sqrt(sparams.hbar/(sparams.me*originOmega));
    CMEsHOs = kron(bcoeffs,bcoeffs)*(CMEs_lib_sub/A)*kron(bcoeffs,bcoeffs)';
    CMEsItin = kron(acoeffs,acoeffs)*CMEsHOs*kron(acoeffs,acoeffs)';

    % Get the basis states that you want to use based on number of electrons in
    % the system
    H2ndQ = buildSecondQuantizationHam(sparams, CMEsItin);

    [manyBody_evecs, manyBody_ens] = eigs(H2ndQ,8,'sa');
end

