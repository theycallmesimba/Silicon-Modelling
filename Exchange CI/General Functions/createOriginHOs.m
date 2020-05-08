function [originHOs, originOmega] = createOriginHOs( sparams, gparams, omega, waitbarFlag )
%CREATENONSHIFTEDHOS Summary of this function goes here
%   Detailed explanation goes here

    if nargin < 4
        waitbarFlag = 1;
    end

    XX = gparams.XX;
    YY = gparams.YY;
    
    % If omega is input as an argument, then use that, otherwise use the
    % result from the potential fitting.
    if nargin < 3
        originOmega = abs(mean(sparams.fittedPotentialParameters(:,1)));
    else
        originOmega = omega;
    end
        
    if waitbarFlag
        h = waitbar(0,'1','Name','Creating Non Shifted 2D HOs...',...
            'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
    end

    originHOs(sparams.nOriginHOs) = twoDimNonShiftHO;
    alpha = sqrt(sparams.me*originOmega/sparams.hbar);
    
    xx = XX(1,:);
    yy = YY(:,1);
    
    kk = 0;
    
    xHermites = zeros(sparams.maxOriginHOsX,length(xx));
    yHermites = zeros(sparams.maxOriginHOsY,length(yy));
    % It is more efficient to store each hermite polynomial as it is
    % constructed since we will use all of them from 0 -> n to build the
    % harmonic orbitals in x and y
    % X Hermite polynomials
    for ii = 0:(sparams.maxOriginHOsX-1)
        switch ii
            case 0
                xHermites(ii+1,:) = getHermiteN(ii, [], [], alpha*xx);
            case 1
                xHermites(ii+1,:) = getHermiteN(ii, [],...
                    [], alpha*xx);
            otherwise
                 xHermites(ii+1,:) = getHermiteN(ii, xHermites(ii,:),...
                    xHermites(ii-1,:), alpha*xx);
        end
    end
    % Y Hermite polynomials
    for ii = 0:(sparams.maxOriginHOsY-1)
        switch ii
            case 0
                yHermites(ii+1,:) = getHermiteN(ii, [], [], alpha*yy);
            case 1
                yHermites(ii+1,:) = getHermiteN(ii, [],...
                    [], alpha*yy');
            otherwise
                 yHermites(ii+1,:) = getHermiteN(ii, yHermites(ii,:),...
                    yHermites(ii-1,:), alpha*yy');
        end
    end  
    
    xHOs = zeros(sparams.maxOriginHOsX,length(xx));
    yHOs = zeros(sparams.maxOriginHOsY,length(yy));
    eX = zeros(1,sparams.maxOriginHOsX);
    eY = zeros(1,sparams.maxOriginHOsX);
    % X Harmonic Orbitals
    for ii = 0:(sparams.maxOriginHOsX-1)
        xHOs(ii+1,:) = 1/sqrt(2^ii*factorial(ii))*(alpha^2/pi)^(1/4)...
            *exp(-alpha^2*xx.^2/2).*xHermites(ii+1,:);
        eX(ii+1) = sparams.hbar*originOmega*(ii + 1/2);
    end
    % Y Harmonic Orbitals
    for ii = 0:(sparams.maxOriginHOsX-1)
        yHOs(ii+1,:) = 1/sqrt(2^ii*factorial(ii))*(alpha^2/pi)^(1/4)...
            *exp(-alpha^2*yy.^2/2).*yHermites(ii+1,:)';
        eY(ii+1) = sparams.hbar*originOmega*(ii + 1/2);
    end
    
    for ii = 0:(sparams.maxOriginHOsX-1)
        for jj = 0:(sparams.maxOriginHOsY-1)
            if waitbarFlag
                %Check for cancel button click
                if getappdata(h,'canceling')
                    break;
                end
            end
            
            WFMG = (yHOs(jj+1,:)'*ones(1,gparams.ngridx)).*(ones(gparams.ngridy,1)*xHOs(ii+1,:));
                        
            kk = kk + 1;
            originHOs(kk).initialize(WFMG,convertMGtoNO(WFMG),ii,jj,eX(ii+1)+eY(jj+1));
            
            if waitbarFlag
                % Update waitbar
                waitbar(kk/(sparams.nOriginHOs), h,...
                    sprintf('Current X index:%d/%d  Current Y index:%d/%d',ii,...
                    sparams.maxOriginHOsX,jj,sparams.maxOriginHOsY));
            end
        end
    end
    if waitbarFlag
        % Close waitbar
        delete(h);
    end
end

function hermiteN = getHermiteN(nn, hermiteNm1, hermiteNm2, xx)
    switch nn
        case 0
            hermiteN = ones(1,length(xx));
        case 1
            hermiteN = 2*xx;
        otherwise
            hermiteN = (2*xx.*hermiteNm1) - (2*(nn-1)*hermiteNm2);
    end
end
