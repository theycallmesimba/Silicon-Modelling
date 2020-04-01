function sparams = calculateStarkShift( sparams, pot, psi, xx, jj, yy )
%CALCULATESTARKSHIFT Summary of this function goes here
%   Detailed explanation goes here

    % Find the electric field
    [~,currEz] = gradient(pot/sparams.ee,sparams.dx,sparams.dz);
    currEz = -currEz;

    % Find ground state of current potential
    groundState = solve1DSingleElectronSE(sparams,1,xx,pot(sparams.twoDEGindZ,:));

    % Find average Ez seen by the ground state and current
    % simulated wavefunction
    sparams.avgEzGround(jj,yy) = getInnerProduct(xx,groundState.',currEz(sparams.twoDEGindZ,:).*groundState.');
    sparams.avgEz(jj,yy) = getInnerProduct(xx,psi,currEz(sparams.twoDEGindZ,:).*psi);
    
    % Find shift in resonance frequency seen by the ground state
    % and the current simulated wavefunction
    sparams.vShiftGround(jj,yy) = sparams.n2*sparams.v0*(abs(sparams.avgEzGround(jj,yy))^2) + sparams.v0;
    sparams.vShift(jj,yy) = sparams.n2*sparams.v0*(abs(sparams.avgEz(jj,yy))^2) + sparams.v0;
end

