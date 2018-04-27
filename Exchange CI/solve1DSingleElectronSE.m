function [ wfs, ens ] = solve1DSingleElectronSE( sparams, X, V )
%SOLVE1DSINGLEELECTRONSE Solves the Schrodinger Equation for a 1D potential
%function using the Finite Difference Method.
%   xx = x coordinates in an array (must be same length as potential and
%   have same grid spacing)
%   potential = potential function values in an array (must be same length
%   as xx)
%   solNum = number of solutions to return starting from the ground state
%   consts = specifies values for calculation (hbar, effective mass,
%   electorn charge, etc)
    
    full1DLap = make1DSELap(sparams,X,V);

    % Determine which eigs option to use to get the correct e-vectors
    [eVectors, ens] = eigs(full1DLap,sparams.nLocalOrbitals,'sa');

    wfs = zeros(length(X),sparams.nLocalOrbitals);
    for kk = 1:sparams.nLocalOrbitals
        for ii = 1:length(X)
            wfs(ii,kk) = eVectors(ii,kk);
        end
    end
end






