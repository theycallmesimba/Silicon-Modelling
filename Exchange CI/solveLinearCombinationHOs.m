function [sparams, HLOHO] = solveLinearCombinationHOs( sparams, X, Y, V )
%GETLINEARCOMBINATIONCOEFFICIENTS Summary of this function goes here
%   Detailed explanation goes here

    % Get the Laplacian
    full2DLap = make2DSELap(sparams,X,Y,V);

    % Now we will rewrite the Hamiltonian in the LOHO basis
    HLOHO = zeros(sparams.nSingleOrbitals);
    SLOHO = zeros(sparams.nSingleOrbitals);
    for ii = 1:(sparams.nSingleOrbitals)
        currWFLeft = sparams.loeLOHOs(ii).wavefunctionMG;
        for jj = 1:(sparams.nSingleOrbitals)
            currWFRight = sparams.loeLOHOs(jj).wavefunctionNO;
            HLOHO(ii,jj) = getInnerProduct2D(currWFLeft,...
                convertNOtoMG(full2DLap*currWFRight,sparams.ngridx,sparams.ngridy),X,Y);
            SLOHO(ii,jj) = getInnerProduct2D(currWFLeft,sparams.loeLOHOs(jj).wavefunctionMG,X,Y);
        end
    end

    % Solve the regular eigenvalue equation
    [vecs,en] = eig(HLOHO,SLOHO);
    [~,ind] = sort(diag(en));
    en = en(ind,ind);
    sparams.acoeffs = vecs(:,ind);
    
    % Transpose matrix to match indexing convention in code
    sparams.acoeffs = sparams.acoeffs.';
    
    % Now build up the single particle wavefunctions using the acoeffs
    sparams.linearCombinationSHOs(sparams.nSingleOrbitals) = twoDimLCHO;
    for ii = 1:sparams.nSingleOrbitals
        tempwf = zeros(sparams.ngridy,sparams.ngridx);
        
        % Now normalize rows of the transformation matrix
        sparams.acoeffs(ii,:) = sparams.acoeffs(ii,:)/norm(sparams.acoeffs(ii,:));

        for jj = 1:sparams.nSingleOrbitals
            tempwf = tempwf + sparams.acoeffs(ii,jj)*sparams.loeLOHOs(jj).wavefunctionMG;
        end
        
        % Normalize the wave function
        tempNorm = sqrt(getInnerProduct2D(tempwf,tempwf,X,Y));
        sparams.LCHOs(ii).wavefunctionMG = tempwf/tempNorm;
        sparams.LCHOs(ii).wavefunctionNO = convertMGtoNO(tempwf)/tempNorm;
        sparams.LCHOs(ii).energy = en(ii,ii);
    end
end
