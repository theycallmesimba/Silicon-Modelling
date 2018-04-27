function sparams = solveLoewdinOrthonormalization( sparams, X, Y )
%SOLVELOEWDINORTHONORMALIZATION Summary of this function goes here
%   Detailed explanation goes here

    % First step is to get the overlap matrix
    smatrixLocalHOs = zeros(sparams.nSingleOrbitals,sparams.nSingleOrbitals);
    for ii = 1:sparams.nSingleOrbitals
        currWFii = sparams.localHOs(ii).wavefunctionMG;
        for jj = 1:sparams.nSingleOrbitals
            currWFjj = sparams.localHOs(jj).wavefunctionMG;
            smatrixLocalHOs(ii,jj) = getInnerProduct(currWFii,currWFjj,X,Y);
        end
    end
    
    % Now get the square root of the matrix
    sparams.smatrixSQRTLocalHOs = sqrtm(smatrixLocalHOs);
    
    sSQRTinv = inv(sparams.smatrixSQRTLocalHOs);
    % Now we want to build up our new basis
    for ii = 1:sparams.nSingleOrbitals
        currWF = zeros(sparams.ngridy,sparams.ngridx);
        for jj = 1:sparams.nSingleOrbitals  
            currWF = currWF + sSQRTinv(ii,jj)*sparams.localHOs(jj).wavefunctionMG;
        end
        sparams.sLocalHOs(ii).wavefunctionMG = currWF;
        sparams.sLocalHOs(ii).wavefunctionNO = convertMGtoNO(currWF);
    end
end

