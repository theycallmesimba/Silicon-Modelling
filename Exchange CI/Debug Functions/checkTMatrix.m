function checkTMatrix( sparams, Tmatrix )
%CHECKTMATRIX Summary of this function goes here
%   Detailed explanation goes here

    isUnitary = 1;

    prod = Tmatrix*Tmatrix';
    diagProd = diag(prod);
    offDiagProd = prod - diag(diagProd);
    
    if any(any(abs(diagProd - 1) >= sparams.normThreshold) == 1)
        isUnitary = 0;
        fprintf(1,'Transformation is not Unitary. Normalization issue.\n');
    end
    if any(any(abs(offDiagProd - 0) >= sparams.normThreshold) == 1)
        isUnitary = 0;
        fprintf(1,'Transformation is not Unitary. Orthogonality issue.\n');
    end
    
    if isUnitary
        fprintf(1,'Transformation is Unitary.\n');
    end
end

