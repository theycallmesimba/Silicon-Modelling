function [corrFid, corrRho, corrPhase] = correctT0rotation(rho)
%CORRECTT0ROTATION Summary of this function goes here
%   Detailed explanation goes here

    dz = [1 0;0 -1];
    rhoT = 1/2*[0; 1; 1; 0]*[0, 1, 1, 0];
    rhoS = 1/2*[0; 1; -1; 0]*[0, 1, -1, 0];
    
    % Correction Hamiltonian (z-rotation on 2nd spin)
    H = kron(eye(8),dz);

    % Calculate amount of rotation into T0 state
    tr = trace(kron(eye(4),rhoT)*rho);
    cphase = acos(1 - 2*abs(tr));
    % Now we don't know if the actual phase is < or > pi so check both
    % cases
    cp1 = cphase;
    cp2 = 2*pi - cphase;
    
    % Correct rho and calculate new fidelity wrt S state
    Ucorr1 = expm(-1i*H*(-cp1)/2);
    Ucorr2 = expm(-1i*H*(-cp2)/2);
    cfid1 = abs(trace(kron(eye(4),rhoS)*Ucorr1*rho*Ucorr1'))^2;
    cfid2 = abs(trace(kron(eye(4),rhoS)*Ucorr2*rho*Ucorr2'))^2;
    
    if cfid1 > cfid2
        corrFid = cfid1;
        corrRho = Ucorr1*rho*Ucorr1';
        corrPhase = cp1;
    else
        corrFid = cfid2;
        corrRho = Ucorr2*rho*Ucorr2';
        corrPhase = cp2;
    end
end

