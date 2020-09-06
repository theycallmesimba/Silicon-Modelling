function [ Tmatrix, HamInNewBasis, ens ] = findTMatrixViaHamiltonian(sparams,...
    gparams, basisToT, nStates )
%FINDTMATRIXVIAHAMILTONIAN Finds a transformation matrix to change between
%   a given basis [basisToT] and the itinerant bass of the natural 
%    Hamiltonian given by VV
%   T: basisToT -> itinerantBasis

% THIS CODE ASSUMES THE BASISTOT IS ORTHOGONAL

    XX = gparams.XX;
    YY = gparams.YY;
    VV = gparams.VV;

    % Get the Laplacian
    full2DLap = make2DSELap(sparams,gparams);

    % Now we will rewrite the Hamiltonian in the desired basis
    HamInNewBasis = zeros(nStates);

    % Find upper triangular elements
    % If a parpool is running, take advantage of it
    if ~isempty(gcp('nocreate'))
        ngridx = gparams.ngridx;
        ngridy = gparams.ngridy;
        parfor jj = 1:nStates
            tempCol = zeros(nStates,1);
            
            currWFRight = basisToT(jj).wavefunctionNO;
            currWFRight = convertNOtoMG(full2DLap*currWFRight, ngridx, ngridy);
            for ii = (jj+1):nStates
                currWFLeft = basisToT(ii).wavefunctionMG;
                tempCol(ii) = getInnerProduct2D(currWFLeft, currWFRight, XX, YY);
            end
            HamInNewBasis(:,jj) = tempCol;
        end
    else
        for jj = 1:nStates
            currWFRight = basisToT(jj).wavefunctionNO;
            currWFRight = convertNOtoMG(full2DLap*currWFRight,gparams.ngridx,gparams.ngridy);
            for ii = (jj+1):nStates
                currWFLeft = basisToT(ii).wavefunctionMG;
                HamInNewBasis(ii,jj) = getInnerProduct2D(currWFLeft,currWFRight,XX,YY);
            %             SHamInNewBasis(ii,jj) = getInnerProduct2D(currWFLeft,basisToT(jj).wavefunctionMG,XX,YY);
            end
        end
    end
    % Find lower triangular elements
    HamInNewBasis = HamInNewBasis + HamInNewBasis';
    % Find diagonal elements
    for ii = 1:nStates
        currWFRight = basisToT(ii).wavefunctionNO;
        currWFRight = convertNOtoMG(full2DLap*currWFRight,gparams.ngridx,gparams.ngridy);
        
        currWFLeft = basisToT(ii).wavefunctionMG;
        HamInNewBasis(ii,ii) = getInnerProduct2D(currWFLeft,currWFRight,XX,YY);
    end
    HamInNewBasis = (HamInNewBasis + HamInNewBasis')/2;
    
    % Solve the regular eigenvalue equation
    [vecs,ens] = eig(HamInNewBasis);
    % Sort e-vectors by increasing e-energy
    [ens,ind] = sort(diag(ens));
    Tmatrix = vecs(:,ind);
    
    % Transpose matrix to match indexing convention in code
    Tmatrix = Tmatrix.';
end 