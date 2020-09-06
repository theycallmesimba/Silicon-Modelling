function inProd = getInnerProduct2D( wfL, wfR, XX, YY )
%GETINNERPRODUCT finds the inner product between two wavefunctions via
%trapezoidal integration
    F = wfL.*wfR;
    dy = YY(2,1) - YY(1,1);
    dx = XX(1,2) - XX(1,1);
    inProd = dy*ftrapz(dx*ftrapz(F',2),1);
end

