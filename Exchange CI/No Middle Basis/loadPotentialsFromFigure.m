function [ sparams, xdata, ydata, zdata ] = loadPotentialsFromFigure( sparams, figName )
%LOADPOTENTIALSFROMFIGURE Summary of this function goes here
%   Detailed explanation goes here
    open(figName);
    h = gcf; % Get current figure
    axesObjs = findobj(h, 'type', 'axes'); % axes handle
    dataObjs = get(axesObjs, 'Children'); % low-level graphics objects in axes
%     objTypes = get(dataObjs, 'Type')
    xdata = get(dataObjs,'XData');
    ydata = get(dataObjs,'YData');
    zdata = get(dataObjs,'ZData');
    close(h);
    
    [sparams, xdata, ydata, zdata] = interpolatePotential(sparams,xdata,ydata,zdata);
end

