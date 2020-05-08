function [ optOmega ] = optimizeOmega( sparams, gparams, omegaGuess, wfsToCompare, tols )
%OPTIMIZEOMEGAORIGINHO This function tries to improve the convergance of
%the many electron energy spectra calcuation by optimize omega chosen to
%form the harmonic orbitals centered at the origin.
%   sparams: simulation parameter file
%   gparams: grid parameter file
%   omegaGuess: starting point for where to search for omega.  A good
%   typical guess is the average harmonic confinement extracted from the
%   potential landscape for each QD.
 
    if nargin < 4
        % We will choose an ideal omega by seeing which choice of omega allows
        % us to best fit the lowest energy itinerant wave functions. It is
        % impractical to fit to a large set of itinerant WFs, so we just choose
        % the first n = number of QDs.
        % Find the first n itinerant wavefunctions
        wfsToCompare = findItinerantBasis(sparams, gparams, sparams.nDots);
    end
    if nargin < 5
        tols = 1E-7;
    end
    
    % Anonymous function so we can pass s(g)params and the itin WFs into
    % the optimizer
    fun = @(w)findWFDifference(w, sparams, gparams, wfsToCompare);

    % Perform the optimization
    % TODO: Figure out a good way to choose tolerance
    options = optimoptions(@fminunc,...%'Display','iter',...
        'Algorithm','quasi-newton','OptimalityTolerance', 1E-7);
    % Sometimes get negative omegas
    optOmega = abs(fminunc(fun,omegaGuess,options));
end

function diffWF = findWFDifference(omegaGuess, sparams, gparams, wfsToCompare)
    % For the current omega guess, construct the basis of harmonic orbitals
    % centered at the origin
    % Switch any negative omegas to positive in the search
    omegaGuess = abs(omegaGuess);
    
    basisNewOmega = createOriginHOs(sparams, gparams, omegaGuess, 0);
    
    % Now calculate the inner product of this basis and the itinerant basis
    normFlag = 0;
    Tmatrix = findTMatrixViaInnerProd(gparams, basisNewOmega, wfsToCompare, normFlag);
    
    % Tmatrix contains all <itin_i|origin_j> inner products.  If we have
    % chosen a good omega, then SUM(|<itin_i|origin_J>|^2) will be close to
    % 1.  We want to maximize this value as it means the basis of origin
    % harmonic orbitals approximates the itin orbitals well.
    minCondition = abs(1-diag(Tmatrix*Tmatrix'));
    % Now average this difference for all of the itin wavefunctions we are
    % considering.
    [~,nWFs] = size(wfsToCompare);
    diffWF = sum(minCondition)/nWFs;
end





