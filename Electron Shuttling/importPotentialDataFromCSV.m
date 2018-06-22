function [pot2D, XXq, ZZq] = importPotentialDataFromCSV(fname)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    data = dlmread(fname);
    [rows,cols] = size(data);
    zdata = data(2:rows,2:cols);
    xdata = data(1,2:cols);
    ydata = data(2:rows,1);
    
    % The data may not be uniform in sampling, so we need to fix that for
    % the fourier transforms in the main code for speed up.
    % Find next highest power of 2 to the length of xx
    desiredGridX = 2^(nextpow2(length(xdata)));
    desiredGridZ = round(2.5*length(ydata));
                    
    % Make linearly spaced grid of points to interpolate the
    % potential at
    xxq = linspace(min(xdata),max(xdata),desiredGridX);
    zzq = linspace(min(ydata),max(ydata),desiredGridZ);

    [XX,ZZ] = meshgrid(xdata,ydata);
    [XXq,ZZq] = meshgrid(xxq,zzq);
    
    pot2D = interp2(XX,ZZ,zdata,XXq,ZZq);
end