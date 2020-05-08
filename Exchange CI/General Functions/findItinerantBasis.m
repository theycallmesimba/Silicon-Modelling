function [itinBasis, ens] = findItinerantBasis( sparams, gparams, nStates )
%FINDITINERANTBASIS Summary of this function goes here
%   Detailed explanation goes here
    
    [wfs,ens] = solve2DSingleElectronSE(sparams, gparams, nStates);
    [ens,perm] = sort(diag(ens));
    wfs = wfs(:,:,perm);
    
    itinBasis = struct;
    for ii = 1:nStates
        itinBasis(ii).wavefunctionMG = squeeze(wfs(:,:,ii));
        itinBasis(ii).wavefunctionNO = convertMGtoNO(squeeze(wfs(:,:,ii)));
    end
end

