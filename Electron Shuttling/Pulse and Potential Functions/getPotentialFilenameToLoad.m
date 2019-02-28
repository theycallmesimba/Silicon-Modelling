function fname = getPotentialFilenameToLoad( sparams, recursionIndex, currFVec, fileFormatType )
%GETPOTENTIALFILENAMETOLOAD Summary of this function goes here
%   Detailed explanation goes here

    if recursionIndex ~= sparams.numOfGates
        currFname = getPotentialFilenameToLoad(sparams, recursionIndex + 1, currFVec, fileFormatType);
        
        fname = [sprintf('V%d_%0.3f_',recursionIndex,...
            sparams.voltagesToLoad{recursionIndex}(currFVec(recursionIndex))),...
            currFname];
    else
        if strcmp(fileFormatType,'csv') 
            fname = sprintf('V%d_%0.3f.csv',recursionIndex,...
                sparams.voltagesToLoad{recursionIndex}(currFVec(recursionIndex)));
        elseif strcmp(fileFormatType,'nextnano')
            fname = sprintf('V%d_%0.3f/output/potential',recursionIndex,...
                sparams.voltagesToLoad{recursionIndex}(currFVec(recursionIndex)));
        end
    end    
end

