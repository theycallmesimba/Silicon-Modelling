function index = getClosestArrayIndex(value, array)
%GETCLOSESTARRAYELEMENTINDEX Summary of this function goes here
%   Detailed explanation goes here
    [~,index] = min(abs(array - value));
end

