function subCMELib = getSubsetOfCMEs( sparams, fullCMELib )
%GETSUBSETOFCMES This function gets a subset of a larger CME library matrix
%which can be useful for debugging convergence.
%   Detailed explanation goes here
    
    nSub = sparams.nOriginHOs;
    
    [nFull, ~] = size(fullCMELib);
    nFull = sqrt(nFull);
    
    % If the number of HOs for the sub matrix matches the number of HOs
    % used to calculate the full matrix, then we can just return the full
    % matrix.
    if nSub == nFull
        subCMELib = fullCMELib;
        return;
    end
        
    subFromFullInd = zeros(nSub^2,1);

    % m first then n
    indexLookup = zeros(nFull,2);
    kk = 1;
    for ii = 1:sqrt(nFull)
        for jj = 1:sqrt(nFull)
            indexLookup(kk,1) = ii;
            indexLookup(kk,2) = jj;
            kk = kk+1;
        end
    end
    
    subInd = 1;
    for rowIndex = 1:nFull^2
        beta = floor((rowIndex-1)/nFull)+1;
        alpha = mod(rowIndex-1,nFull)+1;
        
        nalpha = indexLookup(alpha,1);
        malpha = indexLookup(alpha,2);
        
        nbeta = indexLookup(beta,1);
        mbeta = indexLookup(beta,2);
        
        if nbeta <= sparams.maxOriginHOsX && mbeta <= sparams.maxOriginHOsY &&...
                nalpha <= sparams.maxOriginHOsX && malpha <= sparams.maxOriginHOsY
            subFromFullInd(subInd) = rowIndex;
            subInd = subInd + 1;
        end
    end

    subCMELib = sparse(fullCMELib(subFromFullInd,subFromFullInd));
end

