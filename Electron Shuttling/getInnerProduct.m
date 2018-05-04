function inProd = getInnerProduct( xx, wfL, wfR )
%GETINNERPRODUCT Summary of this function goes here
%   Detailed explanation goes here

    inProd = trapz(xx,conj(wfL).*wfR);
end

