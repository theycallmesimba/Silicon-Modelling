function H = updateEffectiveHamiltonianDetuning(H,depsL,depsR)
%UPDATEEFFECTIVEHAMILTONIANDETUNING Summary of this function goes here
%   Detailed explanation goes here
    [rows,~] = size(H);
    
    % Left dot first
    for ii = 1:(rows/2)
        H(ii,ii) = H(ii,ii) + depsL;
    end
    % Right dot now
    for ii = (rows/2+1):rows
        H(ii,ii) = H(ii,ii) + depsR;
    end
end

