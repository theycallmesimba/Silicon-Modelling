function [ wfs, ens ] = solve1DSingleElectronSE( sparams, nSols, X, V )
%SOLVE1DSINGLEELECTRONSE Solves the Schrodinger Equation for a 1D potential
%function using the Finite Difference Method.
%   X = x coordinates in an array (must be same length as potential and
%   have same grid spacing)
%   V = potential function values in an array (must be same length
%   as xx)
%   solNum = number of solutions to return starting from the ground state
%   consts = specifies values for calculation (hbar, effective mass,
%   electorn charge, etc)
%     [~,sparams.twoDEGindZ] = min(abs(zz - (-0.5*1E-9)));
%     for ii = 1:length(sparams.potentials)
%         sparams.potentials(ii).pot2DEG = sparams.potentials(ii).pot2D(sparams.twoDEGindZ,:);
%     end 
    full1DLap = make1DSELap(sparams,X,V);

    % Determine which eigs option to use to get the correct e-vectors
    [eVectors, ens] = eigs(full1DLap,nSols,'sa');

    wfs = zeros(length(X),nSols);
    for ii = 1:nSols
        wfs(:,ii) = eVectors(:,ii)/sqrt(getInnerProduct(X,eVectors(:,ii),eVectors(:,ii)));
    end
end




