function [sparams, HLOHO, SLOHO] = solveLinearCombinationHOs( sparams, X, Y, V )
%GETLINEARCOMBINATIONCOEFFICIENTS Summary of this function goes here
%   Detailed explanation goes here

    % Get the Laplacian
    full2DLap = make2DSELap(sparams,X,Y,V);

    % Now we will rewrite the Hamiltonian in the LOHO basis
%     HLOHO = zeros(sparams.nSingleOrbitals);
%     SLOHO = zeros(sparams.nSingleOrbitals);
%     for ii = 1:(sparams.nSingleOrbitals)
%         currWFLeft = sparams.localHOs(ii).wavefunctionMG;
%         for jj = 1:(sparams.nSingleOrbitals)
%             currWFRight = sparams.localHOs(jj).wavefunctionNO;
%             HLOHO(ii,jj) = getInnerProduct(currWFLeft,...
%                 convertNOtoMG(full2DLap*currWFRight,sparams.ngridx,sparams.ngridy),X,Y);
%             SLOHO(ii,jj) = getInnerProduct(currWFLeft,sparams.localHOs(jj).wavefunctionMG,X,Y);
%         end
%     end
    HLOHO = zeros(sparams.nSingleOrbitals);
    SLOHO = zeros(sparams.nSingleOrbitals);
    for ii = 1:(sparams.nSingleOrbitals)
        currWFLeft = sparams.sLocalHOs(ii).wavefunctionMG;
        for jj = 1:(sparams.nSingleOrbitals)
            currWFRight = sparams.sLocalHOs(jj).wavefunctionNO;
            HLOHO(ii,jj) = getInnerProduct(currWFLeft,...
                convertNOtoMG(full2DLap*currWFRight,sparams.ngridx,sparams.ngridy),X,Y);
            SLOHO(ii,jj) = getInnerProduct(currWFLeft,sparams.sLocalHOs(jj).wavefunctionMG,X,Y);
        end
    end

    % Solve the regular eigenvalue equation
    [sparams.acoeffs,en] = eig(HLOHO,SLOHO);

    % Order output from smallest to largest
%     [~,perm] = sort(diag(en));
%     en = en(perm,perm);
    sparams.acoeffs = sparams.acoeffs.';
    
    % Now build up the single particle wavefunctions using the acoeffs
    sparams.linearCombinationSHOs(sparams.nSingleOrbitals) = twoDimLCHO;
    for ii = 1:sparams.nSingleOrbitals
        tempwf = zeros(sparams.ngridy,sparams.ngridx);
        
        % Now normalize rows of the transformation matrix
        sparams.acoeffs(ii,:) = sparams.acoeffs(ii,:)/norm(sparams.acoeffs(ii,:));
%         sparams.acoeffs(:,ii) = sparams.acoeffs(:,ii)/norm(sparams.acoeffs(:,ii));


%         for jj = 1:sparams.nSingleOrbitals
%             tempwf = tempwf + sparams.acoeffs(ii,jj)*sparams.localHOs(jj).wavefunctionMG;
%         end
        for jj = 1:sparams.nSingleOrbitals
            tempwf = tempwf + sparams.acoeffs(ii,jj)*sparams.sLocalHOs(jj).wavefunctionMG;
        end
        
        % Normalize the wave function
        tempNorm = sqrt(getInnerProduct(tempwf,tempwf,X,Y));
        sparams.linearCombinationSHOs(ii).wavefunctionMG = tempwf/tempNorm;
        sparams.linearCombinationSHOs(ii).wavefunctionNO = convertMGtoNO(tempwf)/tempNorm;
        sparams.linearCombinationSHOs(ii).energy = en(ii,ii);
    end
end
