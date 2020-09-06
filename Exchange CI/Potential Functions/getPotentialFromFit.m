function [ pots ] = getPotentialFromFit( fit, xx )
%GETPOTENTIALFROMFIT Summary of this function goes here
%   Detailed explanation goes here
    vals = coeffvalues(fit);
    omega = vals(1);
    V = vals(2);
    a = vals(3);
    
    for ii = 1:length(xx)
        pots(ii) = V + (1/4)*omega^2*(xx(ii) - a)^2;
    end
end
