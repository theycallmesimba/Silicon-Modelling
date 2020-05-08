function [factorialLookup, nchoosekLookup, gammaLookup] = createLookupTables( sparams )
%CREATELOOKUPTABLES Summary of this function goes here
%   Detailed explanation goes here

    maxFactorial = max([sparams.maxOriginHOsX,sparams.maxOriginHOsY]);
    factorialLookup = zeros(1,maxFactorial+1); % include 0
    for ii = 0:maxFactorial
        factorialLookup(ii+1) = prod(1:ii);
    end
    
    nchoosekLookup = zeros(maxFactorial+1,maxFactorial+1); % n, k
    for ii = 0:maxFactorial
        for jj = 0:ii
            nchoosekLookup(ii+1,jj+1) = factorialLookup(ii+1)/...
                (factorialLookup(jj+1)*factorialLookup(ii-jj+1));
        end
    end
    
    % The max value that will be inputted into gamma will be when 2p is
    % maximal
    maxGamma = 2*(sparams.maxOriginHOsX + sparams.maxOriginHOsY);
    gammaLookup = zeros(1,maxGamma+1);
    for ii = 0:maxGamma
        gammaLookup(ii+1) = gamma(ii + 0.5);
    end
end

