function inProd = getInnerProduct( wfL, wfR, X, Y )
%GETINNERPRODUCT Summary of this function goes here
%   Detailed explanation goes here

    inProd = trapz(Y(:,1),trapz(X(1,:),wfL.*wfR,2));
end

