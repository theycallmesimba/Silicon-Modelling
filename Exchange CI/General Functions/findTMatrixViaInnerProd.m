function Tmatrix = findTMatrixViaInnerProd( gparams, basis1, basis2, normFlag )
%FINDTMATRIXVIAINNERPROD Finds a transformation matrix to change between
%   a given basis and another basis by calculating the explicit inner
%   products to form the transformations.
%   T: basis1 -> basis2

    if nargin < 4
        normFlag = 1;
    end

    XX = gparams.XX;
    YY = gparams.YY;

    [~,size1] = size(basis1);
    [~,size2] = size(basis2);

    Tmatrix = zeros(size2,size1);
    
    for jj = 1:size1
        state1 = basis1(jj).wavefunctionMG;
        for ii = 1:size2
            state2 = basis2(ii).wavefunctionMG;
            
            % T matrix with T_ij = <alpha_j|r_i>
            Tmatrix(ii,jj) = getInnerProduct2D(state1,state2,XX,YY);
        end
        if flag == 1
            break;
        end
    end
    
    if normFlag
        % Normalize each row 
        for ii = 1:size2
            Tmatrix(ii,:) = Tmatrix(ii,:)/norm(Tmatrix(ii,:));
        end
    end
end

