function h_fig = analyzeRobustness(xVec, yVec, biasDeltaVec, tunDeltaVec, exchangeMat)
%ANALYZEROBUSTNESS Summary of this function goes here
%   Detailed explanation goes here

    exchangeMat = squeeze(exchangeMat);
    if length(size(exchangeMat)) ~= 4
        error('Exchange matrix needs 4 swept variables.');
    end
    
    nXPts = length(biasDeltaVec);
    nYPts = length(tunDeltaVec);

    % Figure out which variable is the x axis
    switch inputname(1)
        case 'Vtun'
            xVecStr = '$V_{\rm tun}$ [Ry*]';
        case 'dotSep'
            xVecStr = 'Dot separation [$a_B$*]';
        case 'lxi'
            xVecStr = 'Dot width [$a_B$*]';
        case 'eccentricity'
            xVecStr = '$l_y/l_x$';
        otherwise
            error('Could not find matching inputname for xVec variable');
    end
    % Figure out which variable is the y axis
    switch inputname(2)
        case 'Vtun'
            yVecStr = '$V_{\rm tun}$ [Ry*]';
        case 'dotSep'
            yVecStr = 'Dot separation [$a_B$*]';
        case 'lxi'
            yVecStr = 'Dot width [$a_B$]';
        case 'eccentricity'
            yVecStr = '$l_y/l_x$';
        otherwise
            error('Could not find matching inputname for yVec variable');
    end

    

    h_fig = figure('Color','white');
    cm = viridis();
    colormap(cm);
    hold on;
    for ii = 1:length(xVec)
        % Build current x limits
        xTempVec = linspace(0,1/length(xVec),nXPts) + (ii-1)/length(xVec);
        for jj = 1:length(yVec)
            % Build current y limits
            yTempVec = linspace(0,1/length(yVec),nYPts) + (jj-1)/length(yVec);

            % Get current grid
            [XX,YY] = meshgrid(xTempVec,yTempVec);
            currExchange = squeeze(exchangeMat(ii,jj,:,:));
%             s = surf(XX,YY,currExchange);
            refExchange = currExchange(ceil(nXPts/2),ceil(nYPts/2));
            deltaJ = abs(currExchange - refExchange);
%             deltaJ = currExchange - refExchange;
            perDiffJ = deltaJ./refExchange;
            s = surf(XX,YY,perDiffJ);
            set(s,'EdgeColor','white','EdgeAlpha',0.5,'FaceColor','interp');
        end
    end
    view(2);
    % Draw overarching thick white lines
    for ii = 1:length(xVec)-1
        line([ii/length(xVec),ii/length(xVec)],[0,1],[1E3,1E3],'Linewidth',5,'Color','white');
    end
    for jj = 1:length(yVec)-1
        line([0,1,],[jj/length(yVec),jj/length(yVec)],[1E3,1E3],'Linewidth',5,'Color','white');
    end
    % Build tick labels
    xTickArray = [];
    xTickLabelArray = {};
    for ii = 1:length(xVec)
        xTickArray(ii) = (ii-1)/length(xVec) + 1/length(xVec)/2;
        xTickLabelArray{ii} = sprintf('%.2f',xVec(ii)); 
    end
    for jj = 1:length(yVec)
        yTickArray(jj) = (jj-1)/length(yVec) + 1/length(yVec)/2;
        yTickLabelArray{jj} = sprintf('%.2f',yVec(jj)); 
    end

    xticks(xTickArray); yticks(yTickArray);
    xticklabels(xTickLabelArray); yticklabels(yTickLabelArray);
    set(gca,'TickLabelInterpreter','latex','Fontsize',14);
    xlabel(xVecStr,'Interpreter','latex','Fontsize',20);
    ylabel(yVecStr,'Interpreter','latex','Fontsize',20);
    
    cb = colorbar();
    cb.TickLabelInterpreter = 'latex';
    cb.FontSize = 14;
    ylabel(cb,'$|\%J|$ change','Interpreter','latex','Fontsize',20);
%     caxis([0,0.1]);
    caxis([0,0.2]);
%     set(gca,'ColorScale','log');
end

