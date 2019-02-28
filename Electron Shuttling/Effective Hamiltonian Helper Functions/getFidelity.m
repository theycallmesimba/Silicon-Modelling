function fidelity = getFidelity(sparams, rho, H)
%GETFIDELITY Summary of this function goes here
%   Detailed explanation goes here
    if ~sparams.includeSecondSpin
        [kets,ens] = eig(H);
        [~, ind] = sort(diag(ens));
        kets = kets(:,ind);
        kets = setKetSameGlobalPhase(kets);

        [rows,~] = size(H);
        currPsiIdeal = zeros(rows,1);
        for gg = 1:length(sparams.stateIndices)
            currPsiIdeal = currPsiIdeal + 1/sqrt(length(sparams.stateIndices))*kets(:,gg);
        end
        currRhoIdeal = currPsiIdeal*currPsiIdeal';

        fidelity = abs(trace(currRhoIdeal*rho))^2;
    else
        singlet = 1/sqrt(2)*[0,1,-1,0]';
        rhoSinglet = singlet*singlet';
        rhoToCheck = kron(eye(4),rhoSinglet);

        fidelity = abs(trace(rhoToCheck*rho))^2;
    end
end

