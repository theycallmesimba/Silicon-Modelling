function [ optOmega ] = optimizeOmega( sparams, gparams, omegaGuess,...
    wfsToCompare )
%OPTIMIZEOMEGAORIGINHO This function tries to improve the convergance of
%the many electron energy spectra calcuation by optimize omega chosen to
%form the harmonic orbitals centered at the origin.
%   sparams: simulation parameter file
%   gparams: grid parameter file
%   omegaGuess: starting point for where to search for omega.  A good
%   typical guess is the average harmonic confinement extracted from the
%   potential landscape for each QD.
 
    % If the grid size is very large and fine... Then this optimization can
    % take a long time. It is useful to create a new grid instead more
    % coarsely spaced as this is just to find an optimal omega which should
    % be the same even with the coarser grid.
    % Get dx and dy for gparams in nm
    if strcmp(sparams.unitsType,'Rydberg')
        [~,dx] = convertRyToSI(sparams, [], gparams.xx(2) - gparams.xx(1));
        [~,dy] = convertRyToSI(sparams, [], gparams.xx(2) - gparams.xx(1));
    end
    gparams_opt = gparams;
    % Check if x-axis has lots of points and a grid spacing < 1 nm
    % If so, then make a coarser x grid
    if gparams_opt.ngridx > 300 && dx < 1E-9
        gparams_opt.ngridx = 300;
        xx_opt = linspace(gparams.xx(1),gparams.xx(end),gparams_opt.ngridx);
        gparams_opt.xx = xx_opt;
    end
    % Check if y-axis has lots of points and a grid spacing < 1 nm
    if gparams_opt.ngridy > 300 && dy < 1E-9
        gparams_opt.ngridy = 300;
        yy_opt = linspace(gparams.yy(1),gparams.yy(end),gparams_opt.ngridx);
        gparams_opt.yy = yy_opt;
    end
    [XX_opt, YY_opt] = meshgrid(gparams_opt.xx,gparams_opt.yy);
    gparams_opt.XX = XX_opt;
    gparams_opt.YY = YY_opt;
    gparams_opt.VV = interp2(gparams.XX,gparams.YY,gparams.VV,...
        gparams_opt.XX,gparams_opt.YY,'spline');

    if nargin < 4
        % We will choose an ideal omega by seeing which choice of omega allows
        % us to best fit the lowest energy itinerant wave functions. It is
        % impractical to fit to a large set of itinerant WFs, so we just choose
        % the first n = number of QDs.
        % Find the first n itinerant wavefunctions
%         wfsToCompare = findItinerantBasis(sparams, gparams_opt, sparams.nDots);
        wfsToCompare = findItinerantBasis(sparams, gparams_opt, 3*sparams.nDots);
    end
    if nargin < 5
        tols = 1E-6;
    end
    
    % Anonymous function so we can pass s(g)params and the itin WFs into
    % the optimizer
    fun = @(w)findWFDifference(w, sparams, gparams_opt, wfsToCompare);

    % Perform the optimization
    % If a parpool is running, take advantage of it
    if ~isempty(gcp('nocreate'))
        options = optimoptions(@fminunc,...%'Display','iter',...
            'Algorithm','quasi-newton','OptimalityTolerance', tols,...
            'UseParallel',true);
    else
        options = optimoptions(@fminunc,...%'Display','iter',...
            'Algorithm','quasi-newton','OptimalityTolerance', tols);
    end
    logOmegaGuess = log10(omegaGuess);
    optOmega = 10^fminunc(fun,logOmegaGuess,options);
end

function diffWF = findWFDifference(logOmegaGuess, sparams, gparams_opt, wfsToCompare)
    % For the current omega guess, construct the basis of harmonic orbitals
    % centered at the origin
    omegaGuess = 10^logOmegaGuess;
    
    basisNewOmega = createOriginHOs(sparams, gparams_opt, omegaGuess, 0);
    
    % Now calculate the inner product of this basis and the itinerant basis
    normFlag = 0;
    Tmatrix = findTMatrixViaInnerProd(gparams_opt, basisNewOmega, wfsToCompare, normFlag);
    
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





