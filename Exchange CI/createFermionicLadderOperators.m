function sparams = createFermionicLadderOperators(sparams)
%CREATEFERMIONICLADDEROPERATORS Summary of this function goes here
%   Detailed explanation goes here
%***CURRENTLY INCORRECT***
    nStates = 2^(2*sparams.nSingleOrbitals); % x2 for both spins
    
    cOp = cell(1,nStates);
    cdagOp = cell(1,nStates);
    
    % Regular ladder operators with only one state possible
    c = [0 1;0 0];
    cdag = [0 0;1 0];
    dZ = [1 0;0 -1]; % Used to make anti-symmetry of fermionic operators
    dI = eye(2); % Used to expand subspace to include all the orbitals
    
    % Basis set we use is as follows...
    % Odd indices create a spin-up electron at orbital (ii+1)/2
    % Even indices create a spin-down electron at orbital ii/2
    for ii = 1:nStates
        tempcOp = sparse(c);
        tempcdagOp = sparse(cdag);
        for jj = 1:nStates
            if ii == jj
                continue;
            elseif ii < jj
                tempcOp = kron(tempcOp,dZ);
                tempcdagOp = kron(tempcdagOp,dZ);
            elseif ii > jj
                tempcOp = kron(dI,tempcOp);
                tempcdagOp = kron(dI,tempcdagOp);
            end
        end
        cOp{ii} = tempcOp;
        cdagOp{ii} = tempcdagOp;
    end
    
    sparams.cOp = cOp;
    sparams.cdagOp = cdagOp;
end

