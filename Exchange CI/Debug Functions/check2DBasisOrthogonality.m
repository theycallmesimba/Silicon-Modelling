function check2DBasisOrthogonality(sparams, gridparams, basisToCheck, nStates )
%CHECK2DBASISSETORTHOGONALITY Summary of this function goes here
%   Detailed explanation goes here
    foundNonOrtho = 0;
        
    XX = gridparams.XX;
    YY = gridparams.YY; 
    
    for ii = 1:nStates
        for jj = 1:nStates
            if ii == jj
                continue
            end
            currNorm = getInnerProduct2D(basisToCheck(ii).wavefunctionMG,...
                basisToCheck(jj).wavefunctionMG,XX,YY);
            if abs(currNorm - 0) >= sparams.normThreshold
                foundNonOrtho = 1;
                fprintf(1,'WARNING: States <%d|%d> are not orthogonal: %g\n',ii,jj,currNorm);
            end   
        end
    end
    
    if ~foundNonOrtho
        fprintf(1,'All states are orthogonal.\n');
    end
end

