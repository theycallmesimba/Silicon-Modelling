function potential = squeezeFast( numOfInterpGates, potential )
%REMOVESINGLETONDIMENSIONS Summary of this function goes here
%   Detailed explanation goes here

    % 1D slice of the potential
    if length(size(potential)) == numOfInterpGates + 1
        potential = permute(potential,[numOfInterpGates+1,1:numOfInterpGates]);
    % 2D slice of the potential
    elseif length(size(potential)) == numOfInterpGates + 2
        potential = permute(potential,[numOfInterpGates+1,numOfInterpGates+2,1:numOfInterpGates]);    
    end
end

