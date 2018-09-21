function fname = getPotentialFilenameToLoad( sparams, recursionIndex, currFVec )
%GETPOTENTIALFILENAMETOLOAD Summary of this function goes here
%   Detailed explanation goes here

    if recursionIndex ~= sparams.numOfGates
        currFname = getPotentialFilenameToLoad(sparams, recursionIndex + 1, currFVec);
        
        fname = [sprintf('V%d_%0.3f_',recursionIndex,...
            sparams.voltagesToLoad{recursionIndex}(currFVec(recursionIndex))),...
            currFname];
    else
        fname = sprintf('V%d_%0.3f.csv',recursionIndex,...
            sparams.voltagesToLoad{recursionIndex}(currFVec(recursionIndex)));
    end    
end

