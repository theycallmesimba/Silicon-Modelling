function h_fig = analyzeBiasSweeps( sparams, biasVec, yVec, exchangeMat, normJFlag, convertUnitsFlag )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    
    % Default to just showing the raw J data in effective Ry and aB units
    if nargin < 5
        normJFlag = 0;
    end
    if nargin < 6
        convertUnitsFlag = 0;
    end

    h_fig = figure('Color','white');
    
    % Convert bias, yvec, and exchange units
    xVecStr = '$V_{bias}$ [Ry*]';
    ZZ = exchangeMat;
    zzVecStr = '$J$ [Ry*]';
    if convertUnitsFlag
        biasVec = convertRyToSI(sparams, biasVec, [])/sparams.ee/1E-3; % meV
        xVecStr = '$V_{bias}$ [meV]';
        ZZ = convertRyToSI(sparams, squeeze(exchangeMat), [])/sparams.ee/1E-6; % mueV
        zzVecStr = '$J$ [$\mu$eV]';
    end
    
    switch inputname(3)
        case 'Vtun'
            ZZ = ZZ';
            yVecStr = '$V_{\rm tun}$ [Ry*]';
            if convertUnitsFlag
                yVec = convertRyToSI(sparams, yVec, [])/sparams.ee/1E-3; % meV
                yVecStr = '$V_{\rm tun}$ [meV]';
            end
        case 'dotSep'
            yVecStr = 'Dot separation [$a_B^*$]';
            if convertUnitsFlag
                [~, yVec] = convertRyToSI(sparams, [], yVec);
                yVec = yVec/1E-9; % nm
                yVecStr = 'Dot separation [nm]';
            end
        case 'lxi'
            yVecStr = 'Dot width [$a_B^*$]';
            if convertUnitsFlag
                [~, yVec] = convertRyToSI(sparams, [], yVec);
                yVec = yVec/1E-9; % nm
                yVecStr = 'Dot width [nm]';
            end
        case 'eccentricity'
            yVec = yVec;
            yVecStr = '$l_y/l_x$';
        otherwise
            inputname(3)
            error('Could not find matching inputname for yVec variable');
            return;
    end
    
    % Normalize J if we need to/can
    % Calculate the over/under rotation from charge noise assuming it is
    % slow wrt exchange pulse (i.e. constant) and that the exchange pulse
    % is a constant square pulse
    if normJFlag
        if ~any(biasVec == 0)
            fprintf(1,'There is no Vbias = 0 to normalize to. Using minimum value as reference.\n');
            J_Vbias0 = repmat(min(ZZ,[],2),1,length(biasVec));
        else
            J_Vbias0 = repmat(ZZ(:,biasVec == 0),1,length(biasVec));
        end
        diffJ = (ZZ - J_Vbias0)./J_Vbias0;
        % Calculate second derivative of \DeltaJ w.r.t bias
        d1 = gradient(diffJ)./gradient(biasVec);
        d2 = gradient(d1)./gradient(biasVec);
%         d3 = gradient(d2)./gradient(biasVec);
%         ZZ = d3;
%         ZZ = d2;
        ZZ = d1;
%         ZZ = diffJ
        
%         bias0Ind = find(biasVec == 0);
%         biasSelectionInd = [1:(bias0Ind-1) (bias0Ind+1):length(biasVec)]; 
%         biasVec = biasVec(biasSelectionInd);
%         ZZ = ZZ(:,biasSelectionInd);
%         
%         ZZ(ZZ <= 0) = NaN;
        
        zzVecStr = '$\partial^2\%J/\partial V_{\rm bias}^2$';
        if convertUnitsFlag
            zzVecStr = '$\partial^2\%J/\partial V_{\rm bias}^2$';
        end
    end
    
    % Form the surface plot
    [XX,YY] = meshgrid(biasVec, yVec);
    s = surf(XX,YY,ZZ,'FaceColor','interp');
    view(2);
    set(s,'EdgeAlpha',0.5,'EdgeColor','w');
    cm = viridis();
    colormap(cm);
    cb = colorbar();
    cb.TickLabelInterpreter = 'latex';
    cb.FontSize = 14;
    ylabel(cb,zzVecStr,'Interpreter','latex','Fontsize',20);
    set(gca,'Fontsize',14,'TickLabelInterpreter','latex');
    xlabel(xVecStr,'Fontsize',20,'Interpreter','latex');
    ylabel(yVecStr,'Fontsize',20,'Interpreter','latex');
    xlim([min(min(XX)),max(max(XX))]);
    ylim([min(min(YY)),max(max(YY))]);
    
%     set(gca,'ColorScale','log');
%     caxis([0,0.1]);
end

