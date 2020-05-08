function [H2ndQ, T, Hc] = buildSecondQuantizationHam( sparams, CMEsItin, debug )
%BUILDSECONDQUANTIZATIONHAM This function builds the second quantization
%Hamiltonian used for solving the many-body Schrodinger equation.
%   sparams: the sparams sctructure which contains a lot of important
%   information regarding the simulation and the coulomb matrix elements
%   required to build the Hamiltonian.
%   CMEsItin: the full set of Coulomb matrix elements in the itinerant
%   basis of the system.

    % Build the spin-orbital basis vectors
    [SOvectors, SOvectors_map] = buildSOBasisVectors(sparams, debug);

    % Now get the current number of states used in the second quantization
    % calculation after truncating the spin sub-space
    [nSOStates,~] = size(SOvectors);
    
    T = sparse(zeros(nSOStates));
    Hc = sparse(zeros(nSOStates));
    
    % The second quantization hamiltonian has two terms: T + H_c.  T
    % describes all the single particle energies in the hamiltonian while
    % H_c describes all the direct and exchange electron interactions.
    
    % Let's first build the single-particle energy operator as it's the
    % easiest.
    % All elements lie only on the diagonal and each element is simply the
    % sum of single electron energies comprising the state
    % sum_j(\eps_j c_j^\dag c_j)
    % First get all of the energies
    itinOrbitalEnergies = sparams.LCHOEnergies;
    
    for ii = 1:nSOStates
        currElecOrbitals = SOvectors(ii,1:sparams.numElectrons);
        T(ii,ii) = sum(itinOrbitalEnergies(currElecOrbitals));
    end
    
%     tic;
    % Now that the T matrix is assembled, let's turn to the H_c term.
%     parfor nn = 1:nSOStates
%         for mm = 1:nSOStates
%             Hc(nn,mm) = hcHelper(sparams, nn, mm, CMEsItin,...
%                 SOvectors, SOvectors_map);
%         end
%     end
%     toc;
    
%     tic;
    % Now that the T matrix is assembled, let's turn to the H_c term.
    for nn = 1:nSOStates
        for mm = 1:nSOStates
            Hc(nn,mm) = hcHelper(sparams, nn, mm, CMEsItin,...
                SOvectors, SOvectors_map);
        end
    end
%     toc;
    
    % Build the full Hamiltonian (and correct numerical errors by forcing
    % it to be symmetric)
    H2ndQ = T + Hc;
    H2ndQ = (H2ndQ + H2ndQ')/2;
end

function hcElem = hcHelper(sparams, nn, mm, CMEsItin, SOvectors, SOvectors_map)

    nElec = sparams.numElectrons;
    
    [nSO,~] = size(SOvectors_map); 

    hcElem = 0;
    
    % Loop over i > j
    %****************%
    for jj = 1:(nSO-1)
        % Initialize the bra state <ij|
        braStateOrig = SOvectors(nn,:);
        
        % jth anhilation operator
        braStateOrig = anhilationOperatorHelper(sparams, braStateOrig, jj, SOvectors_map);
        
        % Used for the final two checks
        jjOrbital = SOvectors_map(jj,1);
        jjSpin = SOvectors_map(jj,2);
        
        for ii = (jj+1):nSO
            % Refresh bra state for new loop
            braState = braStateOrig;
%             % Initialize the bra state <ij|
%             braState = SOvectors(nn,:);
            
            % Apply anhilation operators to BRA state
            % ith anhilation operator
            braState = anhilationOperatorHelper(sparams, braState, ii, SOvectors_map);
            
%             % jth anhilation operator
%             braState = anhilationOperatorHelper(sparams, braState, jj, SOvectors_map);

            % Now that we have the modified bra state, check that both
            % anhilation operators acted on it.  If they both did not, then
            % they can be commuted so the one that did not acts on the
            % vacuum state giving c|0> = 0.
%             n_anhil_applied = length(find(braState == -1));
            n_anhil_applied = sum(braState == -1);
            if n_anhil_applied/2 ~= 2 % Divide by 2 because we modify both the orbital and spin part of braState
                continue;
            end
            
            % Now calculate the phase from applying these anhilation
            % operators
            braPhase = calculatePhaseHelper(braState);
            
            % Used for the final two checks
            iiOrbital = SOvectors_map(ii,1);
            iiSpin = SOvectors_map(ii,2);
            
            % Loop over k < l
            for kk = 1:(nSO-1)
                % Initialize the ket state |kl>
                ketStateOrig = SOvectors(mm,:);
                
                % kth anhilation operator
                ketStateOrig = anhilationOperatorHelper(sparams, ketStateOrig, kk, SOvectors_map);
                
                % Used for the final two checks
                kkOrbital = SOvectors_map(kk,1);
                kkSpin = SOvectors_map(kk,2);
                
                for ll = (kk+1):nSO
                    % Refresh ket state for new loop
                    ketState = ketStateOrig;
%                     % Initialize the ket state |kl>
%                     ketState = SOvectors(mm,:);

                    % Apply anhilation operators to KET state
                    % lth anhilation operator
                    ketState = anhilationOperatorHelper(sparams, ketState, ll, SOvectors_map);
                    
%                     % kth anhilation operator
%                     ketState = anhilationOperatorHelper(sparams, ketState, kk, SOvectors_map);

                    % Now that we have the modified ket state, check that both
                    % anhilation operators acted on it.  If they both did not, then
                    % they can be commuted so the one that did not acts on the
                    % vacuum state giving c|0> = 0.
%                     n_anhil_applied = length(find(ketState == -1));
                    n_anhil_applied = sum(ketState == -1);
                    if n_anhil_applied/2 ~= 2 % Divide by 2 because we modify both the orbital and spin part of ketState
                        continue;
                    end

                    % Now calculate the phase from applying these anhilation
                    % operators
                    ketPhase = calculatePhaseHelper(ketState);
                    
                    % Now that we've done all the modification of the bra
                    % and kets.. Check that all the remaining electrons are
                    % the same for the bra and ket.  Otherwise, the inner
                    % product is 0.
                    braStateMod = braState(braState ~= -1);
                    ketStateMod = ketState(ketState ~= -1);
                    
                    % SANITY CHECK: there should be nElec - 2 electrons remaining.
                    if length(braStateMod)/2 ~= (nElec - 2) ||...
                            length(braStateMod) ~= length(ketStateMod)
                        fprintf(1,'ERROR: incorrect number of electrons left after application of bra/ket annhilation operators.\n');
                        return;
                    end
                    
                    % Now is the check for the remaining electrons
                    % matching.  Because of our ordering convention, after
                    % we take out the anhilated states, the arrays should
                    % be exactly equal if they are the same.  No need to
                    % worry about permutated vectors.
                    if ~isequal(braStateMod, ketStateMod)
                        continue;
                    end
                    
                    % Used for the final two checks
                    llOrbital = SOvectors_map(ll,1);
                    llSpin = SOvectors_map(ll,2);
                    
                    % Check that the spins for state pairs (i,l) and (j,k) 
                    % match
                    if iiSpin == llSpin && jjSpin == kkSpin
                        rowInd = (jjOrbital - 1)*sparams.nItinerantOrbitals + iiOrbital;
                        colInd = (kkOrbital - 1)*sparams.nItinerantOrbitals + llOrbital;
                        currCME = CMEsItin(rowInd,colInd);
                        hcElem = hcElem + currCME*braPhase*ketPhase;
                    end
                    % Check that the spins for state pairs (i,k) and (j,l) 
                    % match
                    if iiSpin == kkSpin && jjSpin == llSpin
                        rowInd = (jjOrbital - 1)*sparams.nItinerantOrbitals + iiOrbital;
                        colInd = (llOrbital - 1)*sparams.nItinerantOrbitals + kkOrbital;
                        currCME = CMEsItin(rowInd,colInd);
                        hcElem = hcElem - currCME*braPhase*ketPhase;
                    end
                end
            end
        end
    end
end

% This function simply takes an input ket and applies the anhilation
% operator corresponding to the n = [ind] spin-orbit state
function state = anhilationOperatorHelper(sparams, state, ind, SOvectors_map)
    
    anhil_ind = SOvectors_map(ind,:);
    nElec = sparams.numElectrons;
    
    for temp = 1:nElec
        % Check if the ith anhilation operator destroys any of the
        % states in the ket
        if state(temp) == anhil_ind(1) && state(temp+nElec) == anhil_ind(2)
            % Edit the orbital and spin state so we know it was destroyed
            state(temp) = -1; state(temp+nElec) = -1;
            % This will only happen once so break out of the loop
            break;
        end
    end
end

% This function takes an inputted ket that has been modified by annhilation
% operators and calculates the phase term caused by swapping of the
% fermionic operators during the annhilation applications
function phase = calculatePhaseHelper(state)

    % First find all the -1 indicies in the state vector as these
    % correspond to anhilated electrons
    anhilated_indices = find(state == -1);
    
    % Truncate the second half of these indices which correspond to spin
    % and not orbital.
    anhilated_indices = anhilated_indices(1:(end/2));
    
    % Subtracting 1 will convert the indicies into how many swaps we had to
    % do to apply each anhilation operator
    required_swaps = anhilated_indices - 1;
    
    % Now calculate the phase
    phase = (-1)^sum(required_swaps);
end




