function [ Tmatrix, HamInNewBasis, ens ] = findTMatrixViaHamiltonian(sparams,...
    gparams, basisToT, nStates )
%FINDTMATRIXVIAHAMILTONIAN Finds a transformation matrix to change between
%   a given basis [basisToT] and the itinerant bass of the natural 
%    Hamiltonian given by VV
%   T: basisToT -> itinerantBasis
    
% TODO: Move the convertNOtoMG in the jj loop to the ii loop to reduces
% redundant calls
    XX = gparams.XX;
    YY = gparams.YY;
    VV = gparams.VV;

    % Get the Laplacian
    full2DLap = make2DSELap(sparams,XX,YY,VV);

    % Now we will rewrite the Hamiltonian in the LoeHO basis
    HamInNewBasis = zeros(nStates);
    SHamInNewBasis = zeros(nStates);
%     for ii = 1:nStates
%         currWFLeft = basisToT(ii).wavefunctionMG;
%         for jj = 1:nStates
%             currWFRight = basisToT(jj).wavefunctionNO;
%             HamInNewBasis(ii,jj) = getInnerProduct2D(currWFLeft,...
%                 convertNOtoMG(full2DLap*currWFRight,gparams.ngridx,gparams.ngridy),XX,YY);
%             SHamInNewBasis(ii,jj) = getInnerProduct2D(currWFLeft,basisToT(jj).wavefunctionMG,XX,YY);
%         end
%     end
    for jj = 1:nStates
        currWFRight = basisToT(jj).wavefunctionNO;
        currWFRight = convertNOtoMG(full2DLap*currWFRight,gparams.ngridx,gparams.ngridy);
        for ii = 1:nStates
            currWFLeft = basisToT(ii).wavefunctionMG;
            HamInNewBasis(ii,jj) = getInnerProduct2D(currWFLeft,currWFRight,XX,YY);
            SHamInNewBasis(ii,jj) = getInnerProduct2D(currWFLeft,basisToT(jj).wavefunctionMG,XX,YY);
        end
    end

    HamInNewBasis = (HamInNewBasis + HamInNewBasis')/2;
    
    % Solve the regular eigenvalue equation
    [vecs,ens] = eig(HamInNewBasis,SHamInNewBasis);
    % Sort e-vectors by increasing e-energy
    [ens,ind] = sort(diag(ens));
    Tmatrix = vecs(:,ind);
    
    % Transpose matrix to match indexing convention in code
    Tmatrix = Tmatrix.';
end

