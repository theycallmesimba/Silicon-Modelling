function [ sparams, Xq, Yq, Vq ] = interpolatePotential( sparams, X, Y, V )
%INTERPOLATEPOTENTIAL Summary of this function goes here
%   Detailed explanation goes here

    % Check if we are in units of nm or m
    if any(X >= 1)
        X = X*1E-9; % Convert from nm to m
        Y = Y*1E-9; % Convert from nm to m
    end
    V = V*sparams.ee; % Convert to J

    % Resample the coordinates to be powers of 2
    sparams.ngridx = 2^(nextpow2(length(X)));
    sparams.ngridy = 2^(nextpow2(length(Y)));

    [Xq, Yq]= meshgrid(linspace(min(X),max(X),sparams.ngridx),linspace(min(Y),max(Y),sparams.ngridy));

    Vq = -interp2(X,Y,V,Xq,Yq);
end

