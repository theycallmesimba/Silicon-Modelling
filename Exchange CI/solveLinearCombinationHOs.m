function sparams = solveLinearCombinationHOs( sparams, X, Y, V )
%GETLINEARCOMBINATIONCOEFFICIENTS Summary of this function goes here
%   Detailed explanation goes here

    % Get the Laplacian
    full2DLap = make2DSELap(sparams,X,Y,V);

%     % Now we will rewrite the Hamiltonian in the Loe LOHO basis
%     HLOHO = zeros(sparams.nSingleOrbitals);
%     for ii = 1:(sparams.nSingleOrbitals)
%         currWFLeft = sparams.sLocalHOs(ii).wavefunctionMG;
%         for jj = 1:(sparams.nSingleOrbitals)
%             currWFRight = sparams.sLocalHOs(jj).wavefunctionNO;
%             HLOHO(ii,jj) = getInnerProduct(currWFLeft,...
%                 convertNOtoMG(full2DLap*currWFRight,sparams.ngridx,sparams.ngridy),X,Y);
%         end
%     end
    % Now we will rewrite the Hamiltonian in the non shifted HO basis
%     HLOHO = zeros(sparams.nNonShiftedHOs);
%     for ii = 1:(sparams.nNonShiftedHOs)
%         currWFLeft = sparams.nonShiftedHOs(ii).wavefunctionMG;
%         for jj = 1:(sparams.nNonShiftedHOs)
%             currWFRight = sparams.nonShiftedHOs(jj).wavefunctionNO;
%             HLOHO(ii,jj) = getInnerProduct(currWFLeft,...
%                 convertNOtoMG(full2DLap*currWFRight,sparams.ngridx,sparams.ngridy),X,Y);
%         end
%     end

    % Solve the regular eigenvalue equation
    [acoeffs,en] = eig(HLOHO);

    % Order output from smallest to largest
    [~,perm] = sort(diag(en));
    en = en(perm,perm);
    sparams.acoeffs = acoeffs(:,perm).';
    
    % Now build up the single particle wavefunctions using the acoeffs
    sparams.linearCombinationSHOs(sparams.nSingleOrbitals) = twoDimLCHO;
    for ii = 1:sparams.nSingleOrbitals
        tempwf = zeros(sparams.ngridy,sparams.ngridx);
        
        % Now normalize rows of the transformation matrix
        sparams.acoeffs(ii,:) = sparams.acoeffs(ii,:)/norm(sparams.acoeffs(ii,:));

        for jj = 1:sparams.nSingleOrbitals
            tempwf = tempwf + sparams.acoeffs(ii,jj)*sparams.sLocalHOs(jj).wavefunctionMG;
        end
        
        
        % Normalize the wave function
        tempNorm = sqrt(getInnerProduct(tempwf,tempwf,X,Y));
        sparams.linearCombinationSHOs(ii).wavefunctionMG = tempwf/tempNorm;
        sparams.linearCombinationSHOs(ii).wavefunctionNO = convertMGtoNO(tempwf)/tempNorm;
        sparams.linearCombinationSHOs(ii).energy = en(ii,ii);
    end
%     
%     % Now build up the single particle wavefunctions using the acoeffs
%     sparams.linearCombinationSHOs(sparams.nNonShiftedHOs) = twoDimLCHO;
%     for ii = 1:sparams.nSingleOrbitals
%         tempwf = zeros(sparams.ngridy,sparams.ngridx);
%         
%         % Now normalize rows of the transformation matrix
%         sparams.acoeffs(ii,:) = sparams.acoeffs(ii,:)/norm(sparams.acoeffs(ii,:));
%         
%         for jj = 1:sparams.nNonShiftedHOs
%             tempwf = tempwf + sparams.acoeffs(ii,jj)*sparams.nonShiftedHOs(jj).wavefunctionMG;
%         end
%         
%         % Normalize the wave function
%         sparams.linearCombinationNSHOs(ii).wavefunctionMG = tempwf;
%         sparams.linearCombinationNSHOs(ii).wavefunctionNO = convertMGtoNO(tempwf);
%         sparams.linearCombinationNSHOs(ii).energy = en(ii,ii);
%     end
end

