function sparams = fitDotLocationsToHarmonicWells( sparams, X, Y, V )
%FITDOTLOCATIONSTOHARMONICWELLS Summary of this function goes here
%   Detailed explanation goes here

    func = strcat(strcat('Vi + (1/2)*',num2str(sparams.me)),'*wi^2*(xx - ai)^2');
    % func = 'Vi + (1/2)*m*wi^2*(xx - ai)^2';
    modelVariables = {'wi','Vi','ai'};
    fmodel = fittype(func,'ind',{'xx'},'coeff',modelVariables);
        
    sparams.fittedPotentialParameters = zeros(sparams.nDots,length(modelVariables)+1); % +1 for both x and y location
    sparams.fittedXPotentials = zeros(sparams.nDots,sparams.ngridx);
    sparams.fittedYPotentials = zeros(sparams.nDots,sparams.ngridy);

    % For each dot location, we want to fit a HO well
    for ii = 1:sparams.nDots
        % Get the index locations for the meshgrid corresponding to the dot's
        % center location
        [~,xTemp] = min(abs(X(1,:) - sparams.dotLocations(ii,1)));
        [~,yTemp] = min(abs(Y(:,1) - sparams.dotLocations(ii,2)));

        % Along the slice of x or y where the dot is located, find all the
        % tunnel barrier peaks (note: can be empty arrays)
        [xTuns, xTunInd] = findpeaks(V(yTemp,:));
        [yTuns, yTunInd] = findpeaks(V(:,xTemp));

        % Find the tunnel barriers that are closest to the current dot of
        % interest
        tTemp = [];
        if ~isempty(xTunInd)
            [~,ind] = min(xTunInd - xTemp);
            tTemp = [tTemp, xTuns(ind)]; 
        end
        if ~isempty(yTunInd)
            [~,ind] = min(yTunInd - yTemp);
            tTemp = [tTemp, yTuns(ind)]; 
        end
        currTunBar = min(tTemp);

        % Now starting, at the center of the dot, we want to go sweep to the
        % sides to find at what index value we reach slightly below the tunnel
        % barrier
        % First X direction
        Vxfit = [];
        Xfit = [];
        for jj = xTemp:-1:1
            if V(yTemp,jj) < currTunBar
                Vxfit = [Vxfit V(yTemp,jj)];
                Xfit = [Xfit X(yTemp,jj)];
            else
                break;
            end
        end
        Vxfit = flip(Vxfit);
        Xfit = flip(Xfit);
        for jj = (xTemp+1):1:sparams.ngridx
            if V(yTemp,jj) < currTunBar
                Vxfit = [Vxfit V(yTemp,jj)];
                Xfit = [Xfit X(yTemp,jj)];
            else
                break;
            end
        end
        % Now Y direction
        Vyfit = [];
        Yfit = [];
        for jj = yTemp:-1:1
            if V(jj,xTemp) < currTunBar
                Vyfit = [Vyfit V(jj,xTemp)];
                Yfit = [Yfit Y(jj,xTemp)];
            else
                break;
            end
        end
        Vyfit = flip(Vyfit);
        Yfit = flip(Yfit);
        for jj = (yTemp+1):1:sparams.ngridy
            if V(jj,xTemp) < currTunBar
                Vyfit = [Vyfit V(jj,xTemp)];
                Yfit = [Yfit Y(jj,xTemp)];
            else
                break;
            end
        end

        % Get initial guesses for omega (if the guess isn't good the fitting
        % function has a hard time getting a correct fit)
        [~,iTemp] = min(abs(Xfit - sparams.dotLocations(ii,1)));
        wxg = sqrt(2/sparams.me*(currTunBar - Vxfit(iTemp)))/abs(Xfit(end) - sparams.dotLocations(ii,1));
        [~,iTemp] = min(abs(Yfit - sparams.dotLocations(ii,2)));
        wyg = sqrt(2/sparams.me*(currTunBar - Vyfit(iTemp)))/abs(Yfit(end) - sparams.dotLocations(ii,2));

        % Do the actual fit
        myfitX = fit(Xfit', Vxfit', fmodel, 'Start',[wxg, min(Vxfit), sparams.dotLocations(ii,1)]);
        myfitY = fit(Yfit', Vyfit', fmodel, 'Start',[wyg, min(Vyfit), sparams.dotLocations(ii,2)]);

    %     figure;
    %     plot(myfitX,Xfit,Vxfit,'bo-');
    %     xlim([Xfit(1),Xfit(end)]);
    %     
    %     figure;
    %     plot(myfitY,Yfit,Vyfit,'ro-');
    %     xlim([Yfit(1),Yfit(end)]);
    %     
        % Extract the fit parameters and construct the fitted potential
        vals = coeffvalues(myfitX);
        wx = vals(1);Vx = vals(2);ax = vals(3);
        vals = coeffvalues(myfitY);
        wy = vals(1);Vy = vals(2);ay = vals(3);
        w = (wx + wy)/2; % Fit to a symmetric harmonic well
        Vmin = (Vx + Vy)/2; % Obviously our fitted well should have a single minimum
        
        % Store the fitted parameters
        sparams.fittedPotentialParameters(ii,:) = [w,Vmin,ax,ay];
        
        sparams.fittedXPotentials(ii,:) = Vmin + (1/2)*sparams.me*w^2*(X(yTemp,:) - ax).^2;
        sparams.fittedYPotentials(ii,:) = (Vmin + (1/2)*sparams.me*w^2*(Y(:,xTemp) - ay).^2)';
    end
    
    % figure;
    % hold on;
    % s = surf(Xq,Yq,Vq);
    % set(s,'edgecolor','none');
    % minVq = min(min(Vq));
    % maxVq = max(max(Vq));
    % dz = maxVq - minVq;
    % zlim([minVq - dz*0.2,maxVq + dz*0.2]);
    % alpha(s,0.4);
    % [ny,nx] = size(Xq);
    % potHO1 = (potHOY1*ones(1,nx) + ones(ny,1)*potHOX1)/2;
    % s1 = surf(Xq,Yq,potHO1);
    % set(s1,'edgecolor','none');
    % potHO2 = (potHOY2*ones(1,nx) + ones(ny,1)*potHOX2)/2;
    % s2 = surf(Xq,Yq,potHO2);
    % set(s2,'edgecolor','none');
    % alpha(s1,0.6);
    % alpha(s2,0.6);
end