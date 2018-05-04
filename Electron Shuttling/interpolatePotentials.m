function [ newPots ] = interpolatePotentials( oldPots, params )
%INTERPOLATEPOTENTIALS Summary of this function goes here
%   Detailed explanation goes here

    finerGrid = 5000; % How much to increase the potential grid
    [nPots,xGrid] = size(oldPots);
    newPots = zeros(nPots*finerGrid,xGrid);
    
    % Generate time index array saying when the potentials should be
    % updated with ones from the simulation
    updatePotTimeInd = round(linspace(1,params.nFrames,nPots));
    desiredPotTimeGrid = round(linspace(1,params.nFrames,nPots*finerGrid));
%     size(updatePotTimeInd)
    % Now, go through each xgrid point and interpolate linearly
    for ii = 1:xGrid
%         fprintf(1,'%d\n',ii)
        newPots(:,ii) = interp1q(updatePotTimeInd',oldPots(:,ii),desiredPotTimeGrid');
    end
end

