function originCMEs_lib = constructOriginCMEsLibGPU( sparams, waitBarFlag )
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
    
    % Initialize the waitbar
    if waitBarFlag
        hWaitBar = waitbar(0,sprintf('CME element row:%06d/%d  col:%06d/%d',1,numOfWFs^2,1,numOfWFs^2),...
            'Name','Calculating the origin CMEs...',...
            'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
    end
    
    for rowIndex = 1:numOfWFs^2
        
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
        
        colIndArr = gpuArray(rowIndex:numOfWFs^2);
        rowIndArr = gpuArray(ones(1,length(colIndArr))*rowIndex);
                
        % Append the found row of CMEs to the full array
        gpuf = @(colInd, rowIndex, numOfWFs) CME_GPUHelper(colInd, rowIndex, numOfWFs);
        originCMEs_lib(rowIndex,colIndArr) = gather(arrayfun(gpuf, colIndArr, rowIndArr, numOfWFs));  
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

function elem = CME_GPUHelper(colIndex, rowIndex, numOfWFs)
    nX_or_Y = sqrt(numOfWFs);
    
    beta = floor((rowIndex-1)/numOfWFs)+1;
    alpha = mod(rowIndex-1,numOfWFs)+1; 

    nalpha = floor((alpha-1)/nX_or_Y)+1;
    malpha = mod(alpha-1,nX_or_Y)+1; 

    nbeta = floor((beta-1)/nX_or_Y)+1;
    mbeta = mod(beta-1,nX_or_Y)+1; 
    
    gamma = floor((colIndex-1)/numOfWFs)+1;
    delta = mod(colIndex-1,numOfWFs)+1; 

    ngamma = floor((gamma-1)/nX_or_Y)+1;
    mgamma = mod(gamma-1,nX_or_Y)+1; 

    ndelta = floor((delta-1)/nX_or_Y)+1;
    mdelta = mod(delta-1,nX_or_Y)+1; 
        
    % Find the CME
    elem = calculateOriginCMEGPU(nalpha, malpha,...
        nbeta, mbeta, ngamma, mgamma, ndelta, mdelta);
end