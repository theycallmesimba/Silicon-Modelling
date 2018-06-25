 function [ sparams, xxq, zzq ] = loadPotentials( sparams )
%LOADPOTENTIALS Function to either generate the potential profiles
%automatically or load them from an external file
    
    % Declare which simulated voltages you wish to load.
    voltages = [0.6,0.693,0.785,0.8];

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
            waitbar(nn/(length(voltages)^sparams.numOfGates),...
                h, sprintf('Loading file: %s',currFName));
        end
        
        % Build the current file name
        currFName = [];
        for ii = 1:sparams.numOfGates
            currFName = [currFName, sprintf('V%d_%0.3f_',ii,voltages(currFileVec(ii)))];
        end
        currFName = [currFName(1:end-1),'.csv'];
        
        % Load the csv file
        data = dlmread([sparams.potDir currFName]);
        [rows,cols] = size(data);
        zdata = data(2:rows,2:cols);

        % If we are on the first potential
        if nn == 1
            xdata = data(1,2:cols);
            ydata = data(2:rows,1);
            
            cutXData = find(xdata >= -200 & xdata <= 200);
            xdata = xdata(cutXData);

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
        end
        % Now interpolate potential and save it along with current voltage
        % values
        sparams.potentials(nn).pot2D = -sparams.ee*interp2(XX,ZZ,zdata(:,cutXData),XXq,ZZq); % Convert to J and invert
        sparams.potentials(nn).gateValues = voltages(currFileVec);
%         sparams.potentials(nn).gateValues = [voltages(currFileVec(3)),voltages(currFileVec(2)),...
%             voltages(currFileVec(1)),voltages(currFileVec(5)),voltages(currFileVec(4))];
        % Since we've saved our potential, increment the counter
        nn = nn + 1;
        
        % Break condition of the loop
        timeToBreak = ones(1,sparams.numOfGates);
        for ii = 1:sparams.numOfGates
            if currFileVec(ii) == length(voltages)
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
            if currFileVec(ii) > length(voltages)
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

 





