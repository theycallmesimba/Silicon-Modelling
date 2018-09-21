function [ wfs, ens ] = solve1DSingleElectronSE( sparams, nSols, X, V )
%SOLVE1DSINGLEELECTRONSE Solves the Schrodinger Equation for a 1D potential
%function using the Finite Difference Method.
%   X = x coordinates in an array (must be same length as potential and
%   have same grid spacing)
%   V = potential function values in an array (must be same length
%   as xx)
%   solNum = number of solutions to return starting from the ground state

    full1DLap = make1DSELap(sparams,X,V);

    % sigma = min(V) means will find the eigenvalues closest to min(V)
    % which is faster than doing the 'sa' option.
    [eVectors, ens] = eigs(full1DLap,nSols,min(V));
    
    % Output is sorted in descending order so let's fix that
    [~, ind] = sort(diag(ens));
    ens = ens(ind,ind);
    eVectors = eVectors(:,ind);
    
    wfs = zeros(length(X),nSols);
    for ii = 1:nSols
        wfs(:,ii) = eVectors(:,ii)/sqrt(getInnerProduct(X,eVectors(:,ii),eVectors(:,ii)));
    end
end




