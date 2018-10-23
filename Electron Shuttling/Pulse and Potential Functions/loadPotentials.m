 function [ sparams, xxq, zzq ] = loadPotentials( sparams )
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
    set(h,'Position',[0.25,0.4,0.3,0.08]);
    
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
        if mod(nn,5) == 0
            waitbar(nn/nPotentialsToLoad, h, sprintf('Loading file: %s',currFName));
        end
        
        % Build the current file name
        currFName = getPotentialFilenameToLoad( sparams, 1, currFileVec );
        
        if nn == 1
            [xdata, ydata, zdata] = loadPotentialFile([sparams.potDir currFName]);
            
            % The data may not be uniform in sampling, so we need to fix that for
            % the fourier transforms in the main code for speed up.
            % Find next highest power of 2 to the length of xx
            desiredGridX = 2^(nextpow2(length(xdata)));
            desiredGridZ = round(2.5*length(ydata));

            % Make linearly spaced grid of points to interpolate the
            % potential at
            xxq = linspace(min(xdata),max(xdata),desiredGridX);
            zzq = linspace(min(ydata),max(ydata),desiredGridZ);

            [XX,ZZ] = meshgrid(xdata,ydata);
            [XXq,ZZq] = meshgrid(xxq,zzq);
        else
            [~, ~, zdata] = loadPotentialFile(sparams, currFName);
        end
        
        % Now interpolate potential and save it along with current voltage
        % values
        sparams.potentials(nn).pot2D = -sparams.ee*interp2(XX,ZZ,zdata,XXq,ZZq); % Convert to J and invert
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
    
    xxq = xxq*1E-9; % Convert to m
    zzq = zzq*1E-9; 
end

 





