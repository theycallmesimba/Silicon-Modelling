function potential = squeezeFast( numOfGates, potential )
%REMOVESINGLETONDIMENSIONS Summary of this function goes here
%   Detailed explanation goes here

    % 1D slice of the potential
    if length(size(potential)) == numOfGates + 1
        potential = permute(potential,[numOfGates+1,1:numOfGates]);
    % 2D slice of the potential
    elseif length(size(potential)) == numOfGates + 2
        potential = permute(potential,[numOfGates+1,numOfGates+2,1:numOfGates]);    
    end
end

