function [ wfs, ens ] = solve2DSingleElectronSE( sparams, gparams, numSols )
%SOLVE2DSINGLEELECTRONSE Solve the 2D Schrodinger equation using the finite
%difference method.
%   X, Y, and V must be inputted as meshgrid coordinates
%   wfs are outputted as meshgrid matrices
%   ens are outputted as a diagonal matrix
    
    full2DLap = make2DSELap(sparams, gparams);
    
    minVV = min(min(gparams.VV));
    [eVectors, ens] = eigs(full2DLap,numSols,minVV);
        
    wfs = zeros(gparams.ngridy,gparams.ngridx,numSols);
    for kk = 1:numSols
        % Make sure the wavefunction is normalized
        eVectors(:,kk) = eVectors(:,kk)/sqrt(getInnerProduct2D(convertNOtoMG(eVectors(:,kk),gparams.ngridx,gparams.ngridy),...
            convertNOtoMG(eVectors(:,kk),gparams.ngridx,gparams.ngridy),gparams.XX,gparams.YY));

        wfs(:,:,kk) = convertNOtoMG(eVectors(:,kk),gparams.ngridx,gparams.ngridy);
    end
end

