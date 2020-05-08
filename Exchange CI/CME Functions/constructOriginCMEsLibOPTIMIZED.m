function originCMEs_lib = constructOriginCMEsLibOPTIMIZED( sparams, waitBarFlag )
%SOLVECMESSAMEORBITAL Summary of this function goes here
%   Detailed explanation goes here
    
    if nargin == 1
        waitBarFlag = 1;
    end

    % Start our parpool if it hasn't been started already
    if isempty(gcp('nocreate'))
        parpool('local',feature('numcores'));
    end

    % We need to loop through all of the possible two electron
    % configurations.  <alpha,beta|v|gamma,delta>
    % This part will be long depending on your basis size.
    numOfWFs = sparams.nOriginHOs;
    originCMEs_lib = sparse(numOfWFs^2,numOfWFs^2);
    
    % Produce lookup table for n and m values depending on what HO index we
    % are on.  This is needed for parallel computation in which we'll get
    % rid of any reference to sparams here since we don't want the 
    % overhead in sending sparams to each worker in our parallel pool.
    indexLookup = zeros(sparams.nOriginHOs,2);
    for ii = 1:sparams.nOriginHOs
        indexLookup(ii,1) = sparams.originHOs(ii).n;
        indexLookup(ii,2) = sparams.originHOs(ii).m;
    end
    
    % One way we will speed up the calculation, is by precalculating all of
    % the factorial terms used to derive the CMEs (this includes nchoosek,
    % gamma, and beta functions as well)
    % Helper function is at the end of this file
    [factorialLookup, nchoosekLookup, gammaLookup, betaLookup] = createLookupTables(sparams);
    
    % Initialize the waitbar
    if waitBarFlag
        hWaitBar = waitbar(0,sprintf('CME element row:%06d/%d  col:%06d/%d',1,numOfWFs^2,1,numOfWFs^2),...
            'Name','Calculating the origin CMEs...',...
            'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
    end
    
    for rowIndex = 1:numOfWFs^2
        beta = floor((rowIndex-1)/numOfWFs)+1;
        alpha = mod(rowIndex-1,numOfWFs)+1; 
        
        nalpha = indexLookup(alpha,1);
        malpha = indexLookup(alpha,2);
        
        nbeta = indexLookup(beta,1);
        mbeta = indexLookup(beta,2);
        
        %***Waitbar stuff***%
        if waitBarFlag
            % Update waitbar
            waitbar(rowIndex/(numOfWFs^2), hWaitBar,...
                sprintf('CME element row:%03d/%d  col:%03d/%d',rowIndex,numOfWFs^2,rowIndex,numOfWFs^2));
            %Check for cancel button click
            if getappdata(hWaitBar,'canceling')
                break;
            end
        end
        
        tempRow = zeros(1,numOfWFs^2);
        parfor colIndex = rowIndex:numOfWFs^2
            gamma = floor((colIndex-1)/numOfWFs)+1;
            delta = mod(colIndex-1,numOfWFs)+1;

            ngamma = indexLookup(gamma,1);
            mgamma = indexLookup(gamma,2);

            ndelta = indexLookup(delta,1);
            mdelta = indexLookup(delta,2);
            
            % Find the CME
%             tempRow(colIndex) = calculateOriginCME(nalpha, malpha,...
%                 nbeta, mbeta, ngamma, mgamma, ndelta, mdelta);
            tempRow(colIndex) = calculateOriginCMEOPTIMIZED(nalpha, malpha,...
                nbeta, mbeta, ngamma, mgamma, ndelta, mdelta,...
                factorialLookup, nchoosekLookup, gammaLookup, betaLookup);
        end
        
        % Append the found row of CMEs to the full array
        originCMEs_lib(rowIndex,:) = tempRow;
    end
    if waitBarFlag
        delete(hWaitBar);
    end
    
    % We only found the upper triangular part of the matrix so find the
    % lower triangular part here
    temp = originCMEs_lib - diag(diag(originCMEs_lib));
    originCMEs_lib = originCMEs_lib + temp';
    
    if strcmp(sparams.unitsType,'SI')
        % TODO: This will be incorrect as the derivation now assumes the
        % Coulomb potential is in Rydberg units = 2/|\vec{r}_1 -
        % \vec{r}_2|.  Before adding capabilities for SI units, please
        % double check the coefficients for Rydberg and SI units to keep
        % consistent.
        % As a final step, we need to multiply our CME by the scaling factor from
        % our Coulomb potential e^2/(4*pi*e_r)
        k = 1/(4*pi*sparams.eps)*sparams.ee*sparams.ee;
    elseif strcmp(sparams.unitsType,'Rydberg')
        % Derivation is done assuming Rydberg Coulomb potential, so scaling
        % should be 1.
        k = 1;
    end
    originCMEs_lib = k*originCMEs_lib;
end

function [factorialLookup, nchoosekLookup, gammaLookup, betaLookup] =...
    createLookupTables( sparams )

    maxFactorial = max([sparams.maxOriginHOsX,sparams.maxOriginHOsY]);
    factorialLookup = zeros(1,maxFactorial+1); % include 0
    for ii = 0:maxFactorial
        factorialLookup(ii+1) = prod(1:ii);
    end
    
    nchoosekLookup = zeros(maxFactorial+1,maxFactorial+1); % n, k
    for ii = 0:maxFactorial
        for jj = 0:ii
            nchoosekLookup(ii+1,jj+1) = factorialLookup(ii+1)/...
                (factorialLookup(jj+1)*factorialLookup(ii-jj+1));
        end
    end
    
    % The max value that will be inputted into gamma will be when 2p is
    % maximal
    maxGamma = 2*(sparams.maxOriginHOsX + sparams.maxOriginHOsY);
    gammaLookup = zeros(1,maxGamma+1);
    for ii = 0:maxGamma
        gammaLookup(ii+1) = gamma(ii + 0.5);
    end
    
    % The max argument for the left term in the beta function (p - (a-1)/2)
    % is when p is maxed and a=0. Hence 2*max(m)
    % The max argument for the right term in the beta function ((a+1)/2) is
    % when a is maxed. Hence 2*max(n)
    maxLeftArg = 2*sparams.maxOriginHOsY;
    maxRightArg = 2*sparams.maxOriginHOsX;
    betaLookup = zeros(maxLeftArg+1,maxRightArg+1); % p-(a-1)/2, a+1/2
    for ii = 0:maxLeftArg
        for jj = 0:maxRightArg
            betaLookup(ii+1,jj+1) = beta(ii + 0.5, jj + 0.5);
        end
    end
end