function [xdata, zdata, pot2Ddata] = loadPotentialFile(sparams, currFName, interpFlag, fileType, interpDirs, trimFlag)
%LOADPOTENTIALFILE Summary of this function goes here
%   Moved this to a separate file because sometimes I need to load a single
%   file and don't want to do a loop over loading a bunch of them

    if nargin < 3
        fileType = 'csv';
        interpDirs = 'xz';
    end

    if strcmp(fileType,'csv')
        [xdata, zdata, pot2Ddata] = helperLoadcsv(currFName);
    elseif strcmp(fileType,'nextnano')
        [xdata, zdata, pot2Ddata] = helperLoadnextnano(currFName);
    end
    
    if interpFlag
        if trimFlag
            % Trim the x axis a little bit
            xCutOffL = getClosestArrayIndex(-110,xdata);
            xCutOffR = getClosestArrayIndex(110,xdata);
            xdata = xdata(xCutOffL:xCutOffR);
            pot2Ddata = pot2Ddata(:,xCutOffL:xCutOffR);
        end
        
        desiredGridX = 2^(nextpow2(length(xdata)));
        desiredGridZ = round(2.5*length(zdata));

        [XX,YY] = meshgrid(xdata,zdata);
                
        % Make linearly spaced grid of points to interpolate the
        % potential at
        if contains(interpDirs,'x')
            xdata = linspace(min(xdata),max(xdata),desiredGridX);
        end
        if contains(interpDirs,'z')
            zdata = linspace(min(zdata),max(zdata),desiredGridZ);
        end
        
        [XXq,YYq] = meshgrid(xdata,zdata);
        
        pot2Ddata = interp2(XX,YY,pot2Ddata,XXq,YYq);
    end
    
    % Convert to nm and J
    xdata = xdata*1E-9;
    zdata = zdata*1E-9;
    pot2Ddata = -sparams.ee*pot2Ddata;
end

% Load the csv file
function [xd, zd, pot2Dd] = helperLoadcsv(currFName)
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
function [xd, zd, pot2Dd] = helperLoadnextnano(currFName)
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
    yIndex = getClosestArrayIndex(0,yd);
    pot2Dd = permute(squeeze(pot3D(:,yIndex,:)),[2,1]);
end





