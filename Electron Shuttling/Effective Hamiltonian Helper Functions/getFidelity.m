function fidelity = getFidelity(sparams, rho, H, rhoToCompare)
%GETFIDELITY Summary of this function goes here
%   Detailed explanation goes here
    
    % If not specified as an argument, default to comparing with the
    % singlet
    if nargin < 4
        ketToCompare = 1/sqrt(2)*[0,1,-1,0]';
        rhoToCompare = ketToCompare*ketToCompare';
    end

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
        rhoToCompare = kron(eye(4),rhoToCompare);
        fidelity = abs(trace(rhoToCompare*rho))^2;
    end
end

