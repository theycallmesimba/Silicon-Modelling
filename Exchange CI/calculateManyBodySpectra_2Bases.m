function [manyBody_evecs, manyBody_ens, CMEs_lib_sub, LCHOEns, optOmega] = calculateManyBodySpectra_2Bases(...
    sparams, gparams, optOmegaFlag, debug, CMEs_lib)
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

    total_sim_time = tic;

    % Check input arguments and assign default values
    if nargin < 5
        CMEs_lib_supplied = 0;
    else
        CMEs_lib_supplied = 1;
    end
    if nargin < 4
        debug = 0;
    end
    if nargin < 3
        optOmegaFlag = 1;
    end

    fprintf(1,'Begining many body energy calculation\n\n');

    %**************************************%
%     omegaGuess = abs(mean(sparams.fittedPotentialParameters(:,1)));
    omegaGuess = sparams.omegaGuess;
    optomega_time = tic;
    if optOmegaFlag
        fprintf(1,'Optimizing origin harmonic orbital omega...\n');
        optOmega = optimizeOmega(sparams,gparams,omegaGuess);
        fprintf(1,'Found an optimal omega of %.2E\n',optOmega);
        toc(optomega_time);
        optomega_time = toc(optomega_time);
        fprintf(1,'Done!\n\n');
    else
        optOmega = omegaGuess;
        optomega_time = toc(optomega_time);
    end
    
    %**************************************%
    fprintf(1,'Finding 2D harmonic orbitals at origin...\n');
    % Create a new basis of HOs centered at the origin of our dots
    [originHOs, originOmega] = createOriginHOs(sparams,gparams,optOmega);

    % Check orthonormality of the states
    if sparams.verbose && debug
        basisToCheck = originHOs;
        nStates = sparams.nOriginHOs;
        plot2DBasis(gparams,basisToCheck,16);
%         check2DBasisOrthogonality(sparams,gparams,basisToCheck,nStates);
%         check2DBasisNormality(sparams,gparams,basisToCheck,nStates);
    end
    fprintf(1,'Done!\n\n');

    %**************************************%
    fprintf(1,'Finding A matrix...\n');
    nStates = sparams.nOriginHOs;
    basisToUse = originHOs;
%     [acoeffs, ~, sparams.LCHOEnergies] = findTMatrixViaHamiltonian(sparams, gparams, basisToUse, nStates);

    aMat_time = tic;
%     itinOrbs = findItinerantBasis(sparams, gparams, sparams.nItinerantOrbitals);
    
%     acoeffs = findTMatrixViaInnerProd(gparams, basisToUse, itinOrbs);
    [ acoeffs, ~, ens ] = findTMatrixViaHamiltonian(sparams,...
        gparams, basisToUse, nStates);
    sparams.LCHOEnergies = ens;

%     % Now find the itinerant basis energies using acoeffs
%     sparams.LCHOEnergies = zeros(1,sparams.nItinerantOrbitals);
%     full2DLap = make2DSELap(sparams,gparams);
%     for ii = 1:sparams.nItinerantOrbitals
%         tempwf = zeros(gparams.ngridy*gparams.ngridx,1);
% 
%         for jj = 1:nStates
%             tempwf = tempwf + acoeffs(ii,jj)*originHOs(jj).wavefunctionNO;
%         end
%         
%         sparams.LCHOEnergies(ii) = getInnerProduct2D(itinOrbs(ii).wavefunctionMG,...
%             convertNOtoMG(full2DLap*tempwf,gparams.ngridx,gparams.ngridy), gparams.XX, gparams.YY);
%     end
    
    % Truncate A in accordance with how many itinerant orbitals we want
    acoeffs = acoeffs(1:sparams.nItinerantOrbitals,:);
    sparams.LCHOEnergies = sparams.LCHOEnergies(1:sparams.nItinerantOrbitals);
    LCHOEns = sparams.LCHOEnergies;
    toc(aMat_time);
    aMat_time = toc(aMat_time);
    
    % Check if the transformation is orthonormal
    if sparams.verbose && debug
        checkTMatrix(sparams, acoeffs);
    end
    % Check the change of basis between the Horigin states and the itinerant
    % basis
    if sparams.verbose && debug
        nStatesToCompare = sparams.nItinerantOrbitals;
        basisToT = originHOs;
        basisToCompare = findItinerantBasis(sparams, gparams, nStatesToCompare);

        checkBasisTransformation(sparams, gparams, basisToT, basisToCompare, acoeffs);
    end
    fprintf(1,'Done!\n\n');
    
    %%%%%%%%%%
%     manyBody_evecs = 0;
%     manyBody_ens = 0;
%     CMEs_lib_sub = CMEs_lib;
%     return;
    %%%%%%%%%%
    
    %**************************************%
    if CMEs_lib_supplied == 1
        % CME library supplied as argument, so no need to load anything
        fprintf(1,'CMEs_lib supplied as an argument!\n');
    else
        fprintf(1,'CMEs_lib not supplied as an argument! Will attempt to load from file...\n');
        % CME library not supplied as an argument, so we need to either
        % load the library if it exists or calculate it.
        if exist(sparams.CMEsLib_fPath, 'file') == 2
            % CME library exists.
            fprintf(1,'Found the CMEs_lib file!\n');
            % Library needs to be loaded
            fprintf(1,'Loading CME library for calculation...  ');
            load(sparams.CMEsLib_fPath);
            assignin('base','CMEs_lib',CMEs_lib);
        else
            % CME library does not exist.
            fprintf(1,'Could not find CMEs_lib file!\n');
            fprintf(1,'Evaluating CMEs for non shifted harmonic orbitals...  \n');
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
    end
    fprintf(1,'Done!\n\n');

    transform_time = tic;
    fprintf(1,'Transforming the CME library to itinerant basis...\n');
    % Get a subset of the CMEs if the given nOrigins parameter in the
    % simparams file is less than what we've solved for in the library
    % already. Useful for checking convergence wrt number of orbitals. 
    CMEs_lib_sub = getSubsetOfCMEs(sparams, CMEs_lib);
    
    % Scale the CMEs to match the origin HOs used to form transformation
    % matrix
    A = sqrt(sparams.hbar/(sparams.me*originOmega));
%     acoeffsTol = 1E-10;
%     acoeffs(acoeffs < acoeffsTol & acoeffs > -acoeffsTol) = 0;
%     acoeffs = sparse(acoeffs);
    kronAcoeffs = kron(acoeffs,acoeffs);
    CMEsItin = (kronAcoeffs*(full(CMEs_lib_sub)/A))*kronAcoeffs';
    toc(transform_time);
    transform_time = toc(transform_time);
    fprintf(1,'Done!\n');
    
    build2nd_time = tic;
    fprintf(1,'Building 2nd quantization Hamiltonian and diagonalizing...\n');
    % Build the 2nd quantization Hamiltonian and then diagonalize it to
    % obtain the egienvectors and eigenenergies of the system.
    debug = 0;
    fprintf(1,'WARNING: DEBUG FOR HERE IS TURNED OFF MANUALLY RIGHT NOW!\n');
    H2ndQ = buildSecondQuantizationHam(sparams, CMEsItin, debug);
    toc(build2nd_time);
    build2nd_time = toc(build2nd_time);
    
    [manyBody_evecs, manyBody_ens] = eigs(H2ndQ,sparams.nOutputtedEnergies,'sa');
    fprintf(1,'Done!\n');
    
    total_sim_time = toc(total_sim_time);
    
    % Display runtime information/etc.
    fprintf(1,'\nTotal simulation time: %.2f sec, %.2f min.\n',...
        total_sim_time,total_sim_time/60);
    fprintf(1,'Time optimizing omega: %.2f sec, %.2f%% of total.\n',...
        optomega_time,optomega_time/total_sim_time*100);
    fprintf(1,'Time finding A matrix: %.2f sec, %.2f%% of total.\n',...
        aMat_time,aMat_time/total_sim_time*100);
    fprintf(1,'Time transforming CMEs: %.2f sec, %.2f%% of total.\n',...
        transform_time,transform_time/total_sim_time*100);
    fprintf(1,'Time building 2nd quantization H: %.2f sec, %.2f%% of total.\n',...
        build2nd_time,build2nd_time/total_sim_time*100);
end

