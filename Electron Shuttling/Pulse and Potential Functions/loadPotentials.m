 function [ sparams, xx, zz ] = loadPotentials( sparams, interpFlag, fileFormatType, interpDirs, trimFlag )
%LOADPOTENTIALS Function to either generate the potential profiles
%automatically or load them from an external file

    nPotentialsToLoad = 1;
    for ii = 1:sparams.numOfGates
        nPotentialsToLoad = nPotentialsToLoad*length(sparams.voltagesToLoad{ii});
    end
        
    % Make the waitbar to show run time
    h = waitbar(0,sprintf('Loading file:'),...
        'Name',sprintf('Loading potentials...'),...
        'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
    set(findall(h,'type','text'),'Interpreter','none');
    set(findall(h),'Units','normalized');
    set(h,'Position',[0.25,0.4,0.5,0.15]);
    
    currFileVec = ones(1,sparams.numOfGates);
    nn = 1;
    
    % Loop will go through every files as given in the voltages array given
    % above (indicating which voltage values to load potentials from) and
    % via the sparams.numOfGates parameter.  Loop will automatically scale
    % and load potentials depending on length of array and number of gates.
    while(1)
        % Check for cancel button click
        if getappdata(h,'canceling')
            flag = 1;
            break;
        end    
        
        % Update waitbar every N loaded potentials
        if mod(nn,3) == 0
            waitbar(nn/nPotentialsToLoad, h, sprintf('Loading file: %s',currFName));
        end
        
        % Build the current file name
        currFName = getPotentialFilenameToLoad( sparams, 1, currFileVec, fileFormatType );
        
        % Load the actual file (interpolation is done during loading here)
        if nn == 1
            [xx, zz, pot2D_XZ] = loadPotentialFile(sparams,[sparams.potDir currFName],...
                interpFlag,fileFormatType,interpDirs,trimFlag);
        else
            [~, ~, pot2D_XZ] = loadPotentialFile(sparams,[sparams.potDir currFName],...
                interpFlag,fileFormatType,interpDirs,trimFlag);
        end
        
        % Convert to J and invert
        sparams.potentials(nn).pot2D = pot2D_XZ;
        currVgVec = [];
        for ii = 1:sparams.numOfGates
            currVgVec = [currVgVec, sparams.voltagesToLoad{ii}(currFileVec(ii))];
        end
        sparams.potentials(nn).gateValues = currVgVec;

        % Since we've saved our potential, increment the counter
        nn = nn + 1;
        
        % Break condition of the loop
        timeToBreak = ones(1,sparams.numOfGates);
        for ii = 1:sparams.numOfGates
            if currFileVec(ii) == length(sparams.voltagesToLoad{ii})
                timeToBreak(ii) = 1;
            else
                timeToBreak(ii) = 0;
            end
        end
        if ~any(timeToBreak == 0)
            break;
        end
        
        % Increment the indexing vector
        for ii = sparams.numOfGates:-1:1
            currFileVec(ii) = currFileVec(ii) + 1;
            if currFileVec(ii) > length(sparams.voltagesToLoad{ii})
                currFileVec(ii) = 1;
            else
                break;
            end
        end
    end
    delete(h);
end

 





