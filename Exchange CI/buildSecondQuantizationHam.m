function sparams = buildSecondQuantizationHam(sparams, spinSubSystems)
%BUILDSECONDQUANTIZATIONHAM Summary of this function goes here
%   Detailed explanation goes here
    nSOstates = 2*sparams.nSingleOrbitals;
    if sparams.numElectrons > nSOstates
        error('Not enough states for the desired number of electrons.\n');
    elseif sparams.numElectrons < 2
        error('Need at least two electrons to calculate J.\n');
    end
    
    % Get all possible state configurations
    stateConfigurations = nchoosek(1:nSOstates,sparams.numElectrons);
    [nStates,~] = size(stateConfigurations);
    
    % Now arrange the states into our format
    sparams.basisVectors = zeros(nStates,2*sparams.numElectrons);
    for ii = 1:nStates
        currentConfig = stateConfigurations(ii,:);
        for jj = 1:sparams.numElectrons
            if mod(currentConfig(jj),2) == 1   
                sparams.basisVectors(ii,jj) = (currentConfig(jj) + 1)/2;
                sparams.basisVectors(ii,jj + sparams.numElectrons) = 0;
            else
                sparams.basisVectors(ii,jj) = currentConfig(jj)/2;
                sparams.basisVectors(ii,jj + sparams.numElectrons) = 1;                
            end
        end
    end
     
    % Now we want to truncate the spin subspace if desired
    % First calculate the possible spin subspaces
    possibleSz = -sparams.numElectrons*0.5:1:sparams.numElectrons*0.5;
    if ischar(spinSubSystems)
        if strcmpi(spinSubSystems,'all')
            % Do nothing as we want the whole spin subspace
        else
            error('Please specify a valid spin subspace system.\n');
        end
    else
        desiredSz = possibleSz(spinSubSystems);
        % Loop through each configuration and calcuate it's Sz.  If it is
        % not one of our desired Sz subspaces, then remove the basis vector.
        for jj = nStates:-1:1
            currSpinConfiguration = sparams.basisVectors(jj,sparams.numElectrons + 1:end);
            currConfigSz = 0;
            for kk = 1:sparams.numElectrons
                if currSpinConfiguration(kk) == 0
                    currConfigSz = currConfigSz - 0.5;
                else
                    currConfigSz = currConfigSz + 0.5;
                end
            end
            if ~any(desiredSz == currConfigSz)
                sparams.basisVectors(jj,:) = [];
            end
        end
    end
    
    [nStates,~] = size(sparams.basisVectors);
    
    % Display all the basis states
    if sparams.verbose
        spinUp = '\x2191';
        spinDown = '\x2193';
        for ii = 1:nStates
            
            basisVecsStr = sprintf('State %d: |',ii);
            for jj = 1:sparams.nSingleOrbitals
                basisVecsStr = [basisVecsStr '0,'];
            end
            basisVecsStr = [basisVecsStr(1:end-1), '>'];
            
            % Used for correctly modifying string when we modify the string
            % and need extra characters
            for jj = sparams.nSingleOrbitals:-1:1
                ind = find(sparams.basisVectors(ii,1:sparams.numElectrons) == jj);
                if length(ind) == 2
                    basisVecsStr = [basisVecsStr(1:(end-2*jj)) spinUp...
                        spinDown basisVecsStr((end-2*jj+2):end)];
                elseif length(ind) == 1
                    spinState = sparams.basisVectors(ii,sparams.numElectrons+ind);
                    if spinState == 0
                        elecSpin = spinDown;
                    else
                        elecSpin = spinUp;
                    end
                    basisVecsStr = [basisVecsStr(1:(end-2*jj)) elecSpin...
                        basisVecsStr((end-2*jj+2):end)];
                end
            end
            fprintf(1,[basisVecsStr, '\n']);
        end
    end
    
    T = sparse(zeros(nStates));
    Hc = sparse(zeros(nStates));
    % The second quantization hamiltonian has two terms: T + H_C.  T
    % describes all the single particle energies in the hamiltonian while
    % H_C describes all the exchange interactions.
    
    % Let's first build the single-particle energy operator as it's the
    % easiest.
    % All elements lie only on the diagonal and each element is simply the
    % sum of single electron energies comprising the state
    % sum_j(\eps_j c_j^\dag c_j)
    % First get all of the energies
    singleElectronOrbitalEnergies = [sparams.LCHOs.energy];
    for ii = 1:nStates
        currElecOrbitals = sparams.basisVectors(ii,1:sparams.numElectrons);
        singleElectronOrbitalEnergies(currElecOrbitals)
        T(ii,ii) = sum(singleElectronOrbitalEnergies(currElecOrbitals));
    end
end

