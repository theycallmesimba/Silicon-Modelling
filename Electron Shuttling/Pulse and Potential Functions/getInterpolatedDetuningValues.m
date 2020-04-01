function [epsL,epsR] = getInterpolatedDetuningValues( qVec, interpolants)
%GETINTERPOLATEDDETUNINGVALUES Summary of this function goes here
%   Detailed explanation goes here
  
    currInterpolant = interpolants{1};
    epsL = currInterpolant({qVec});
    
    currInterpolant = interpolants{2};
    epsR = currInterpolant({qVec});
end
