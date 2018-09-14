function inProd = getInnerProduct2D( wfL, wfR, XX, YY )
%GETINNERPRODUCT finds the inner product between two wavefunctions via
%trapezoidal integration

    inProd = trapz(YY(:,1),trapz(XX(1,:),wfL.*wfR,2));
end

