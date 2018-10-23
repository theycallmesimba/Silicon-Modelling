function [xdata, ydata, zdata] = loadPotentialFile(currFName)
%LOADPOTENTIALFILE Summary of this function goes here
%   Moved this to a separate file because sometimes I need to load a file
%   without doing all the interpotation stuff I do in loadPotentials()

    % Load the csv file
    try
        data = dlmread(currFName);
        [rows,cols] = size(data);
        zdata = data(2:rows,2:cols);

        xdata = data(1,2:cols);
        ydata = data(2:rows,1);
    catch
        error('Failed to load: %s\n',currFName);
    end
end

