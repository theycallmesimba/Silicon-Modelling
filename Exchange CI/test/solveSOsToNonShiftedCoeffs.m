function sparams = solveSOsToNonShiftedCoeffs( sparams, X, Y, V )
%SOLVESOSTONONSHIFTEDCOEFFS Summary of this function goes here
%   Detailed explanation goes here

    % Get the Laplacian
    full2DLap = make2DSELap(sparams,X,Y,V);
    
    % Now we will rewrite the Hamiltonian in the LOHO basis and find the
    % overlap matrix S
    HLOHO = zeros(sparams.nSingleOrbitals);
    SLOHO = zeros(sparams.nSingleOrbitals);
    for ii = 1:(sparams.nSingleOrbitals)
        for jj = 1:(sparams.nSingleOrbitals)
            currWFLeft = sparams.localHOs(ii).wavefunctionNO;
            currWFRight = sparams.localHOs(jj).wavefunctionNO;
            
            HLOHO(ii,jj) = currWFLeft'*full2DLap*currWFRight;
            SLOHO(ii,jj) = currWFLeft'*currWFRight;
        end
    end

    % Solve the generalized eigenvalue problem
    [acoeffs,en] = eig(HLOHO,SLOHO);

    % Order output from smallest to largest
    [~,perm] = sort(diag(en));
    en = en(perm,perm);
    sparams.acoeffs = acoeffs(:,perm).';
    
    % Now build up the single particle wavefunctions using the acoeffs
    sparams.linearCombinationHOs(sparams.nSingleOrbitals) = twoDimLCHO;
    for ii = 1:sparams.nSingleOrbitals
        tempwf = zeros(sparams.ngridy,sparams.ngridx);
        
        for jj = 1:sparams.nSingleOrbitals
            tempwf = tempwf + sparams.acoeffs(ii,jj)*sparams.localHOs(jj).wavefunctionMG;
        end
        
        % Now normalize all the coefficients and wave functions
        tmpNorm = norm(convertMGtoNO(tempwf));
        
        sparams.acoeffs(ii,:) = sparams.acoeffs(ii,:)/tmpNorm;
        
        % Normalize the wave function
        sparams.linearCombinationHOs(ii).wavefunctionMG = tempwf/tmpNorm;
        sparams.linearCombinationHOs(ii).wavefunctionNO = convertMGtoNO(tempwf)/tmpNorm;
        sparams.linearCombinationHOs(ii).energy = en(ii,ii);
    end
end

