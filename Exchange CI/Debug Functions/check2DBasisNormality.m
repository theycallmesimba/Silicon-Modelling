function check2DBasisNormality(sparams, gridparams, basisToCheck, nStates )
%CHECK2DBASISNORMALITY Summary of this function goes here
%   Detailed explanation goes here
    foundNonNorm = 0;
        
    XX = gridparams.XX;
    YY = gridparams.YY; 
    
    for ii = 1:nStates
        currNorm = getInnerProduct2D(basisToCheck(ii).wavefunctionMG,...
            basisToCheck(ii).wavefunctionMG,XX,YY);
        % We will only print out the norm if it is below our acceptable
        % norm tolerance threshold.  Since we are doing things
        % numerically there is bound to be some error in our
        % orthogonality not being explicitly 0
        if abs(currNorm - 1) >= sparams.normThreshold
            foundNonNorm = 1;
            fprintf(1,'WARNING: State |%d> is not normalized: %g\n',ii,currNorm);
        end   
    end
    
    if ~foundNonNorm
        fprintf(1,'All states are normalized.\n');
    end
end

