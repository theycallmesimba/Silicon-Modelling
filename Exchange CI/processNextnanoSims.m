%% SWEEP 2 VARIABLES

pathToOutputFolder = 'C:\Users\bbuonaco\Documents\nextnano\Output\';
pathToSaveFolder = 'C:\Users\bbuonaco\Desktop\Brandon\Exchange Simulations\Sweep_Symscreen_30nm_ecc1.0_ox15_TGate20_0.5nmGrid\';

% fTemplate = 'Vtun_VTUN_%s\\Vbias_BIAS_%s\\AA - DQD_ExFile_VTUN_%s_BIAS_%s';
% fTemplate = 'Vtun_VTUN_%s\\Vbias_BIAS_%s\\AA - DQD_SymScreen_40nm_ecc1.0_ox15_TGate30_0.5nmGrid_VTUN_%s_BIAS_%s';
fTemplate = 'Vtun_VTUN_%s\\AA - DQD_Symscreen_30nm_ecc1.0_ox15_TGate20_0.5nmGrid_VTUN_%s_BIAS_%s';

% sweep1 = 0.116451; % Tunnel gate
sweep1 = 0.03:0.001:0.150;
sweep1SNFlag = 1; % Is sweep 1 in scientific notation?
sweep2 = 0;
% sweep2 = [0,logspace(-7,-2,51)]; % Bias 
% sweep2 = round(sweep2,4,'significant'); % Used for log bias sweeps
sweep2SNFlag = 1; % Is sweep 2 in scientific notation?

% fTemplate = 'DSize_Plg_width_x_%s\\Vtun_VTUN_%s\\AA - DQD_ExFile_Template_Plg_width_x_%s_VTUN_%s';
% sweep1 = [40,50,60,70,80]; % Plunger gate width
% sweep1SNFlag = 0; % Is sweep 1 in scientific notation?
% sweep2 = 0.06:0.0005:0.150; % Tunnel gate 
% sweep2SNFlag = 1; % Is sweep 2 in scientific notation?

% WARNING: CHECK THAT YOU DON"T NEED TO CHANGE THE SAVE STRING

for ii = 1:length(sweep1)
    if sweep1SNFlag
        sweep1Exp = floor(log10(sweep1(ii)));
        sweep1Val = sweep1(ii)/10^sweep1Exp;
        if isnan(sweep1Val)
            sweep1Exp = 0;
            sweep1Val = 0;
        end
        % If number is integer, need to use %.1f formatting for string
        if abs(round(sweep1Val) - sweep1Val) <= 1E-10
            sweep1Str = sprintf('%.1fE%d', sweep1Val, sweep1Exp);
        % Otherwise use %g
        else
            sweep1Str = sprintf('%gE%d', sweep1Val, sweep1Exp);
        end
    else
        sweep1Str = sprintf('%g', sweep1(ii));
        if sweep1(ii) == 0.104550
            sweep1Str = [sweep1Str,'0'];
        end
    end
    
    for jj = 1:length(sweep2)
        if sweep2SNFlag
            sweep2Exp = floor(log10(sweep2(jj)));
            sweep2Val = sweep2(jj)/10^sweep2Exp;
            if isnan(sweep2Val)
                sweep2Exp = 0;
                sweep2Val = 0;
            end
            % If number is integer, need to use %.1f formatting for string
            if abs(round(sweep2Val) - sweep2Val) <= 1E-10
                sweep2Str = sprintf('%.1fE%d', sweep2Val, sweep2Exp);
            % Otherwise use %g
            else
                sweep2Str = sprintf('%gE%d', sweep2Val, sweep2Exp);
            end
        else
            sweep2Str = sprintf('%g', sweep2(jj));
        end
        
        fprintf(1,'Opening sweep file: %d/%d (1) %d/%d (2)\n',...
            ii,length(sweep1),jj,length(sweep2));

       % Build current folder path for sweep conditions
        foldStr = sprintf(fTemplate, sweep1Str, sweep2Str,...
            sweep1Str, sweep2Str); 

        % Build the full file path of interest
        currFName = [pathToOutputFolder, foldStr,...
            '\output\bias_000_000_000_000_000_000_000\potential'];

        % Get file identifier
        [fileID, msg] = fopen([currFName '.coord'],'r');
        
        % Load coordinate data file
        fileData = textscan(fileID,'%f');
        coord = fileData{1};
        fclose(fileID);
        
        % Parse data to get individual x,y,z coords
        coordDiff = diff(coord);
        transitionIndicies = find(coordDiff < 0);

        xd = coord(1:transitionIndicies(1));
        yd = coord(transitionIndicies(1)+1:transitionIndicies(2));
        zd = coord(transitionIndicies(2)+1:end);
        
        % Load potential data and reshape into a 3D matrix
        fileID = fopen([currFName '.dat']);
        fileData = textscan(fileID,'%f');
        pot3D = fileData{1};
        fclose(fileID);
        pot3D = reshape(pot3D,[length(xd),length(yd),length(zd)]);
        
        % Find the 2D potential slice along z where we will assume the 2DEG exists 
        depth_2DEG = -1;
        [~,xy2DEGInd] = min(abs(zd - depth_2DEG));
        pot2DEG = -squeeze(pot3D(:,:,xy2DEGInd))';
        
        savedMatrix = zeros(length(yd)+1,length(xd)+1);
        savedMatrix(1,2:end) = xd';
        savedMatrix(2:end,1) = yd';
        savedMatrix(2:end,2:end) = pot2DEG;
        
        % Now we want to save just the 2D potential to a file
        saveFileTemp = 'Uxy_Vtun_%.6EV_Vbias_%.6EV.txt';
        saveFileName = sprintf(saveFileTemp,sweep1(ii),sweep2(jj));
        
%         pathToSaveFolder = sprintf(pathToSaveFolder_temp);
        dlmwrite([pathToSaveFolder,saveFileName],savedMatrix,'precision','%.13f');
%         writematrix(savedMatrix,[pathToSaveFolder,saveFileName]);
    end
end

%% SWEEP 1 VARIABLE

pathToOutputFolder = 'C:\Users\bbuonaco\Documents\nextnano\Output\';
pathToSaveFolder = 'C:\Users\bbuonaco\Desktop\Brandon\Exchange Simulations\Sweep_Symscreen_80nm_ecc1.0_ox15_TGate20_0.5nmGrid\';

fTemplate = 'Vtun_VTUN_%s\\AA - DQD_Symscreen_80nm_ecc1.0_ox15_TGate20_0.5nmGrid_VTUN_%s';

% sweep1 = .06:0.001:0.150; % Tunnel gate
sweep1 = 0.03:0.001:0.150;
sweep1SNFlag = 1; % Is sweep 1 in scientific notation?
% WARNING: CHECK THAT YOU DON"T NEED TO CHANGE THIS VALUE WHICH IS USED FOR
% SAVING THE DATA.
Vbias = 0;

for ii = 1:length(sweep1)
    if sweep1SNFlag
        sweep1Exp = floor(log10(sweep1(ii)));
        sweep1Val = sweep1(ii)/10^sweep1Exp;
        if isnan(sweep1Val)
            sweep1Exp = 0;
            sweep1Val = 0;
        end
        % If number is integer, need to use %.1f formatting for string
        if abs(round(sweep1Val) - sweep1Val) <= 1E-10
            sweep1Str = sprintf('%.1fE%d', sweep1Val, sweep1Exp);
        % Otherwise use %g
        else
            sweep1Str = sprintf('%gE%d', sweep1Val, sweep1Exp);
        end
    else
        sweep1Str = sprintf('%d', sweep1(ii));
    end
        
    fprintf(1,'Opening sweep file: %d/%d (1) \n',ii,length(sweep1));

   % Build current folder path for sweep conditions
    foldStr = sprintf(fTemplate, sweep1Str, sweep1Str); 

    % Build the full file path of interest
%     currFName = [pathToOutputFolder, foldStr,...
%         '\output\bias_000_000_000_000_000_000_000\potential'];
    currFName = [pathToOutputFolder, foldStr,...
        '\output\bias_000_000_000_000_000_000_000_000\potential'];

    % Get file identifier
    [fileID, msg] = fopen([currFName '.coord'],'r');

    % Load coordinate data file
    fileData = textscan(fileID,'%f');
    coord = fileData{1};
    fclose(fileID);

    % Parse data to get individual x,y,z coords
    coordDiff = diff(coord);
    transitionIndicies = find(coordDiff < 0);

    xd = coord(1:transitionIndicies(1));
    yd = coord(transitionIndicies(1)+1:transitionIndicies(2));
    zd = coord(transitionIndicies(2)+1:end);

    % Load potential data and reshape into a 3D matrix
    fileID = fopen([currFName '.dat']);
    fileData = textscan(fileID,'%f');
    pot3D = fileData{1};
    fclose(fileID);
    pot3D = reshape(pot3D,[length(xd),length(yd),length(zd)]);

    % Find the 2D potential slice along z where we will assume the 2DEG exists 
    depth_2DEG = -1;
    [~,xy2DEGInd] = min(abs(zd - depth_2DEG));
    pot2DEG = -squeeze(pot3D(:,:,xy2DEGInd))';

    savedMatrix = zeros(length(yd)+1,length(xd)+1);
    savedMatrix(1,2:end) = xd';
    savedMatrix(2:end,1) = yd';
    savedMatrix(2:end,2:end) = pot2DEG;

    % Now we want to save just the 2D potential to a file
    saveFileTemp = 'Uxy_Vtun_%.6EV_Vbias_%.6EV.txt';
    saveFileName = sprintf(saveFileTemp,sweep1(ii),Vbias);
    dlmwrite([pathToSaveFolder,saveFileName],savedMatrix,'precision','%.13f');
%         writematrix(savedMatrix,[pathToSaveFolder,saveFileName]);
end


