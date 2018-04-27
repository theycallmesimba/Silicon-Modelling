function [ wfs, ens ] = solve2DSingleElectronSE( sparams, X, Y, V, numSols )
%SOLVE2DSINGLEELECTRONSE Solve the 2D Schrodinger equation using the finite
%difference method.
%   X, Y, and V must be inputted as meshgrid coordinates
%   wfs are outputted as meshgrid matrices
%   ens are outputted as a diagonal matrix
    full2DLap = make2DSELap(sparams,X,Y,V);

    [eVectors, ens] = eigs(full2DLap,numSols,'sa');
    
    [ny,nx] = size(X);
    wfs = zeros(ny,nx,numSols);
    for kk = 1:numSols
        % Make sure the wavefunction is normalized
        eVectors(:,kk) = eVectors(:,kk)/norm(eVectors(:,kk));

        wfs(:,:,kk) = convertNOtoMG(eVectors(:,kk),nx,ny);
    end
end

