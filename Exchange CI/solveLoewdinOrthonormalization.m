function sparams = solveLoewdinOrthonormalization( sparams, X, Y )
%SOLVELOEWDINORTHONORMALIZATION Summary of this function goes here
%   Detailed explanation goes here

    % First step is to get the overlap matrix
    smatrixLOHOs = zeros(sparams.nSingleOrbitals,sparams.nSingleOrbitals);
    for ii = 1:sparams.nSingleOrbitals
        currWFii = sparams.LOHOs(ii).wavefunctionMG;
        for jj = 1:sparams.nSingleOrbitals
            currWFjj = sparams.LOHOs(jj).wavefunctionMG;
            smatrixLOHOs(ii,jj) = getInnerProduct2D(currWFii,currWFjj,X,Y);
        end
    end
    smatrixLOHOs
    % Now get the square root of the matrix
    sSQRTinv = sqrtm(smatrixLOHOs);
    
    sSQRTinv = inv(sSQRTinv);
    % Now we want to build up our new basis
    for ii = 1:sparams.nSingleOrbitals
        currWF = zeros(sparams.ngridy,sparams.ngridx);
        for jj = 1:sparams.nSingleOrbitals  
            currWF = currWF + sSQRTinv(ii,jj)*sparams.LOHOs(jj).wavefunctionMG;
        end
        sparams.loeLOHOs(ii).wavefunctionMG = currWF;
        sparams.loeLOHOs(ii).wavefunctionNO = convertMGtoNO(currWF);
    end
end

