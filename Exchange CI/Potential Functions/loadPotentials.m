function [ xx, pots ] = loadPotentials( fname, consts )
%LOADPOTENTIALS Function to either generate the potential profiles
%automatically or load them from an external file
    
    % Load data from file
    ddir = 'potentials/';
%     data = csvread(strcat(ddir,fname),1,0);
    data = xlsread(strcat(ddir,fname));
    
    % First row of pots should be the position values
    unsampledxx = data(1,:)*1E-9; % Convert to nm
    
    % Rest of the rows are the potentials
    [rows,~] = size(data);
    unsampledpots = data(2:rows,:);
    
    % The data may not be uniform in sampling, so we need to fix that for
    % the fourier transforms in the main code
    % Find next highest power of 2 to the length of xx
    desiredGrid = 2^(nextpow2(length(unsampledxx))+2);
    xlength = max(unsampledxx) - min(unsampledxx);
    % Get the sampling frequency.  Resample assumes that the points before
    % and after your function are 0, so we "pad" our sampling frequency
    % such that the accrued oscillations can be neglected and we still
    % achieve a power of 2 grid spacing.
    nPad = 20;
    dx = xlength/(desiredGrid+(2*nPad));
    [temppots, tempxx] = resample(unsampledpots, unsampledxx, 1/dx,'spline'); 
    pots = temppots(nPad+2:(length(temppots)-nPad-1))*consts.ee; % Convert to Joule's
    xx = tempxx(nPad+2:(length(tempxx)-nPad-1))';
end

