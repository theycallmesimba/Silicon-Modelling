function [xdata, ydata, pot2Ddata] = loadPotentialFileAlt(sparams, currFName,...
    interpFlag, fileType, interpDirs)
%LOADPOTENTIALFILE Summary of this function goes here
%   Moved this to a separate file because sometimes I need to load a single
%   file and don't want to do a loop over loading a bunch of them

    if nargin < 3
        fileType = 'csv';
        interpDirs = 'xy';
    end

    if strcmp(fileType,'csv')
        [xdata, ydata, pot2Ddata] = helperLoadcsv(currFName);
    elseif strcmp(fileType,'nextnano')
        [xdata, ydata, pot2Ddata] = helperLoadnextnano(currFName);
    end
    
    if interpFlag
        
        desiredGridX = 2^(nextpow2(length(xdata)));
        desiredGridY = 2^(nextpow2(length(ydata)));

        [XX,YY] = meshgrid(xdata,ydata);
                
        % Make linearly spaced grid of points to interpolate the
        % potential at
        if contains(interpDirs,'x')
            xdata = linspace(min(xdata),max(xdata),desiredGridX);
        end
        if contains(interpDirs,'y')
            ydata = linspace(min(ydata),max(ydata),desiredGridY);
        end
        
        [XXq,YYq] = meshgrid(xdata,ydata);
        
        pot2Ddata = interp2(XX,YY,pot2Ddata,XXq,YYq);
    end
    
    % Convert to nm and J
    xdata = xdata*1E-9;
    ydata = ydata*1E-9;
    pot2Ddata = -sparams.ee*pot2Ddata;
end

% Load the csv file
function [xd, yd, pot2Dd] = helperLoadcsv(currFName)
    try
        data = dlmread(currFName);
        [rows,cols] = size(data);
        pot2Dd = data(2:rows,2:cols);

        xd = data(1,2:cols);
        zd = data(2:rows,1);
    catch
        error('Failed to load csv file: %s\n',currFName);
    end
end

% Load the nextnano file
function [xd, yd, pot2Dd] = helperLoadnextnano(currFName)
    % First load the coord file to get the correct number of data points
    % for each axis
    coord = dlmread([currFName '.coord']);
    coordDiff = diff(coord);
    transitionIndicies = find(coordDiff < 0);

    xd = coord(1:transitionIndicies(1));
    yd = coord(transitionIndicies(1)+1:transitionIndicies(2));
    zd = coord(transitionIndicies(2)+1:end);

    pot3D = dlmread([currFName '.dat']);
    pot3D = reshape(pot3D,[length(xd),length(yd),length(zd)]);

    % Get the xz potential slice along y = 0 nm
    zIndex = getClosestArrayIndex(-1,zd);
    pot2Dd = squeeze(pot3D(:,:,zIndex))';
end





