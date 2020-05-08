function [corrFid, corrRho] = correctSingletYrotation(rho)
%CORRECTSINGLETXROTATION Summary of this function goes here
%   Detailed explanation goes here

    dy = [0 -1i;1i 0];
    rhoUU_DD = 1/2*[1; 0; 0; 1]*[1, 0, 0, 1];
    rhoS = 1/2*[0; 1; -1; 0]*[0, 1, -1, 0];
    
    % Correction Hamiltonian (x-rotation on 1st spin)
    H = kron(kron(eye(4),dy),eye(2));

    % Calculate amount of rotation into |UU>-|DD> state (which is what |S>
    % rotates into if X is applied to first spin
    tr = trace(kron(eye(4),rhoUU_DD)*rho);
    cphase = asin(sqrt(abs(tr)))*2;
    
    % Now we don't know if the actual angle is < or > pi so check both
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
%         corrPhase = cp1;
        corrRho = Ucorr1*rho*Ucorr1';
    else
        corrFid = cfid2;
%         corrPhase = cp2;
        corrRho = Ucorr2*rho*Ucorr2';
    end
end

