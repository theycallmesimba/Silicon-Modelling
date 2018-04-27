function sparams = createNonShiftedHOs( sparams, X, Y )
%CREATENONSHIFTEDHOS Summary of this function goes here
%   Detailed explanation goes here

    h = waitbar(0,'1','Name','Creating Non Shifted 2D HOs...',...
        'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');

    sparams.nonShiftedHOs((sparams.maxNonShiftedHOsX+1)*(sparams.maxNonShiftedHOsY+1)) = twoDimNonShiftHO;
    omega = abs(mean(sparams.fittedPotentialParameters(:,1)));
    alpha = sparams.me*omega/sparams.hbar;
    
    xx = X(1,:);
    yy = Y(:,1);
    
    kk = 0;
    for ii = linspace(0,sparams.maxNonShiftedHOsX,sparams.maxNonShiftedHOsX+1)
        xWF = 1/sqrt(2^ii*factorial(ii))*(alpha/pi)^(1/4)...
            *exp(-alpha*xx.^2/2).*hermiteH(ii,sqrt(alpha)*xx);
        for jj = linspace(0,sparams.maxNonShiftedHOsY,sparams.maxNonShiftedHOsY+1)
            %Check for cancel button click
            if getappdata(h,'canceling')
                break;
            end
            
            yWF = 1/sqrt(2^jj*factorial(jj))*(alpha/pi)^(1/4)...
                *exp(-alpha*yy.^2/2).*hermiteH(jj,sqrt(alpha)*yy);
            
            WFMG = (yWF*ones(1,sparams.ngridx)).*(ones(sparams.ngridy,1)*xWF);
            WFNO = convertMGtoNO(WFMG);
            normOfWF = norm(WFNO);
            
            kk = kk + 1;
            sparams.nonShiftedHOs(kk).initialize(WFMG/normOfWF,WFNO/normOfWF,ii,jj);
            
            % Update waitbar
            waitbar(kk/((sparams.maxNonShiftedHOsX+1)*(sparams.maxNonShiftedHOsY+1)), h,...
                sprintf('Current X index:%d/%d  Current Y index:%d/%d',ii,...
                sparams.maxNonShiftedHOsX+1,jj,sparams.maxNonShiftedHOsY+1));
        end
    end
    % Close waitbar
    delete(h);
end

