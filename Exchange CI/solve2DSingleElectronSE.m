function [ wfs, ens ] = solve2DSingleElectronSE( sparams, XX, YY, VV, numSols )
%SOLVE2DSINGLEELECTRONSE Solve the 2D Schrodinger equation using the finite
%difference method.
%   X, Y, and V must be inputted as meshgrid coordinates
%   wfs are outputted as meshgrid matrices
%   ens are outputted as a diagonal matrix
    
    full2DLap = make2DSELap(sparams,XX,YY,VV);
    
    minVV = min(min(VV));
    [eVectors, ens] = eigs(full2DLap,numSols,minVV);
        
    [ny,nx] = size(XX);
    wfs = zeros(ny,nx,numSols);
    for kk = 1:numSols
        % Make sure the wavefunction is normalized
        eVectors(:,kk) = eVectors(:,kk)/sqrt(getInnerProduct2D(convertNOtoMG(eVectors(:,kk),nx,ny),...
            convertNOtoMG(eVectors(:,kk),nx,ny),XX,YY));

        wfs(:,:,kk) = convertNOtoMG(eVectors(:,kk),nx,ny);
    end
end

