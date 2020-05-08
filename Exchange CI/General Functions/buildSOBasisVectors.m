function [SOBasisVectors, SOBasisVectors_map] = buildSOBasisVectors( sparams, debug )
%BUILDSOBASISVECTORS This function constructs all of the many body spin 
%orbit basis vectors
%   sparams: details of simulation parameters initialized in simparams.m
%   debug: flag either 0 or 1 to initialize debug mode

    % Check to see if spinsubspaces array is valid or not (S_z indices
    % can't exceed nElectrons+1 or be <= 0.  If invalid, then default to
    % 'all'
    spinSubSpaces = sparams.spinSubspaces;
    if (any(spinSubSpaces > (sparams.numElectrons+1)) ||...
            any(spinSubSpaces <= 0)) && isnumeric(spinSubSpaces)
        fprintf(1,'Invalid S_z subspaces defined. Defaulting to "all".\n');
        spinSubSpaces = 'all';
    elseif ~strcmp(spinSubSpaces,'all') && ~strcmp(spinSubSpaces,'debug') &&...
            ~isnumeric(spinSubSpaces)
        fprintf(1,'Invalid S_z subspaces defined. Defaulting to "all".\n');
        spinSubSpaces = 'all';
    end

    % Get total number of spin-orbit states
    nSOStates = 2*sparams.nItinerantOrbitals;
    if sparams.numElectrons > nSOStates
        error('Not enough states for the desired number of electrons.\n');
    elseif sparams.numElectrons < 2
        error('Need at least two electrons to calculate J.\n');
    end
    
    % Get all possible state configurations and total number (will be cut
    % down a bit later)
    stateConfigurations = nchoosek(1:nSOStates,sparams.numElectrons);
    [n2QStates,~] = size(stateConfigurations);
    
    % Here we create a map between the ith spin-orbit state and the
    % corresponding explicit orbital and spin state.
    SOBasisVectors_map = zeros(nSOStates,2);
    for ii = 1:nSOStates
        SOBasisVectors_map(ii,:) = [floor((ii + 1)/2), mod(ii+1,2)];
    end
    % Order the spin-orbital basis by spin then orbital
    SOBasisVectors_map = sortrows(SOBasisVectors_map,[2,1]);
    
    % Now decode the states found using nchoosek into our format where the
    % first N = sparams.numElectrons indicies correspond to the orbital 
    % state and the last N = sparams.numElectrons indicies correspond to 
    % the spin state [4,2,3,0,0,1] where N = 3 means the first electron is
    % in the 4th orbital state with spin down, the second electron is in 
    % the 2nd orbital state with spin down, and the third electron is in 
    % the 3rd orbital state with spin up.
    SOBasisVectors = zeros(n2QStates,2*sparams.numElectrons);
    for ii = 1:n2QStates
        curConfig = stateConfigurations(ii,:);
        curVec = zeros(1,2*sparams.numElectrons);
        for jj = 1:sparams.numElectrons
            curVec(jj) = SOBasisVectors_map(curConfig(jj),1);
            curVec(jj+sparams.numElectrons) = SOBasisVectors_map(curConfig(jj),2);
        end
        SOBasisVectors(ii,:) = curVec;
    end
    
    % Rewrite each state to the following convention:
    % From left to right, first spin-down and then spin-up.  Within a spin
    % species from left to right, sort by ascending energy level (i.e.
    % orbital index i)
    for ii = 1:n2QStates
        tempVector = SOBasisVectors(ii,:);
        tempVector = reshape(tempVector,[sparams.numElectrons,2]);
        tempVector = sortrows(tempVector,[2,1]);
        SOBasisVectors(ii,:) = reshape(tempVector,[2*sparams.numElectrons,1]);
    end
     
    % Now we want to truncate the spin subspace if desired
    % First calculate the possible spin subspaces
    possibleSz = -sparams.numElectrons*0.5:1:sparams.numElectrons*0.5;
    if ischar(spinSubSpaces)
        % Do nothing as we want the whole spin subspace
        if strcmpi(spinSubSpaces,'all')
        % Debug mode can be used to define basis vectors manually
        elseif strcmpi(spinSubSpaces,'debug')
            SOBasisVectors = [];
            SOBasisVectors(1,:) = [1,2,1,0,0,1];
            SOBasisVectors(2,:) = [2,3,2,0,0,1];
        end
    else
        desiredSz = possibleSz(spinSubSpaces);
        % Loop through each configuration and calcuate Sz.  If it is not
        % one of the desired Sz subspaces, then remove that basis vector.
        % Note that we loop backwards through the states here because we
        % remove rows from the basisVectors matrix which readjusts the
        % indicies of the non-deleted rows.
        for jj = n2QStates:-1:1
            currSpinConfiguration = SOBasisVectors(jj,sparams.numElectrons + 1:end);
            currConfigSz = 0;
            for kk = 1:sparams.numElectrons
                if currSpinConfiguration(kk) == 0
                    currConfigSz = currConfigSz - 0.5;
                else
                    currConfigSz = currConfigSz + 0.5;
                end
            end
            if ~any(desiredSz == currConfigSz)
                SOBasisVectors(jj,:) = [];
            end
        end
    end
    
    % Display all the basis states
    if debug
        printManyElectronBasisVectors(sparams, SOBasisVectors);
    end
end

function printManyElectronBasisVectors( sparams, SOBasisVectors )
%PRINTMANYELECTRONBASISVECTORS Simple function to display the explicit
%forms of the basis vectors for the many electron hamiltonian.

    spinUp = '\x2191';
    spinDown = '\x2193';
    [nStates, ~] = size(SOBasisVectors);
    fprintf(1,'\n');
    for ii = 1:nStates

        basisVecsStr = sprintf('State %d: |',ii);
        for jj = 1:sparams.nItinerantOrbitals
            basisVecsStr = [basisVecsStr '0,'];
        end
        basisVecsStr = [basisVecsStr(1:end-1), '>'];

        % Used for correctly modifying string when we modify the string
        % and need extra characters
        for jj = sparams.nItinerantOrbitals:-1:1
            ind = find(SOBasisVectors(ii,1:sparams.numElectrons) == jj);
            if length(ind) == 2
                basisVecsStr = [basisVecsStr(1:(end-2*jj)) spinUp...
                    spinDown basisVecsStr((end-2*jj+2):end)];
            elseif length(ind) == 1
                spinState = SOBasisVectors(ii,sparams.numElectrons+ind);
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
