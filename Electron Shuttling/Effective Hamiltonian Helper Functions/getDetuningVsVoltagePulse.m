function [epsL, epsR] = getDetuningVsVoltagePulse( sparams, xx, vPulse, showWaitbar )
%GETDETUNINGVSVOLTAGEPULSE Summary of this function goes here
%   Detailed explanation goes here
    [~,nPts] = size(vPulse);

    epsL = zeros(1,nPts);
    epsR = zeros(1,nPts);
    
    if showWaitbar
        h = waitbar(0,sprintf('Current point index: %d/%d',0,nPts),...
            'Name','Finding detuning versus voltage pulse...',...
            'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
    end
    
    for ii = 1:nPts
        % Check for cancel button click
        if showWaitbar
            if getappdata(h,'canceling')
                break;
            end
        end
        if mod(ii,20) == 0 && showWaitbar
            waitbar(ii/nPts, h, sprintf('Current point index: %d/%d',ii,nPts));
        end
        currPotential = sparams.P2DEGInterpolant(getInterpolantArgument(vPulse(:,ii),xx));
        currPotential = squeezeFast(sparams.numOfGates,currPotential)';

        [~, epsL(ii), epsR(ii)] = calculateTunnelCoupling( sparams, xx, currPotential );
    end
    if showWaitbar
        delete(h);
    end
    
    % Account for outliers by doing a constant extrapolation
    startInd = 0;
    foundZero = 0;
    interpPoints = 0;
    for ii = 1:nPts
        if epsL(ii) == 0 && ~foundZero
            startInd = ii;
            foundZero = 1;
        end
        if epsL(ii) ~= 0 && foundZero
            endInd = ii - 1;
            foundZero = 0;
            interpPoints = 1;
        elseif ii == nPts && foundZero
            endInd = ii;
            foundZero = 0;
            interpPoints = 1;
        end
        
        if interpPoints
            interpPoints = 0;
            
            if startInd == 1
                epsL(startInd:endInd) = ones(1,length(startInd:endInd))*epsL(endInd+1);
            elseif endInd == nPts
                epsL(startInd:endInd) = ones(1,length(startInd:endInd))*epsL(startInd-1);
            else
                epsL(startInd:endInd) = interp1([startInd-1,endInd+1],[epsL(startInd-1),epsL(endInd+1)],startInd:endInd);
            end
        end
    end
    
    startInd = 0;
    foundZero = 0;
    interpPoints = 0;
    for ii = 1:nPts
        if epsR(ii) == 0 && ~foundZero
            startInd = ii;
            foundZero = 1;
        end
        if epsR(ii) ~= 0 && foundZero
            endInd = ii - 1;
            foundZero = 0;
            interpPoints = 1;
        elseif ii == nPts && foundZero
            endInd = ii;
            foundZero = 0;
            interpPoints = 1;
        end
        
        if interpPoints
            interpPoints = 0;
            
            if startInd == 1
                epsR(startInd:endInd) = ones(1,length(startInd:endInd))*epsR(endInd+1);
            elseif endInd == nPts
                epsR(startInd:endInd) = ones(1,length(startInd:endInd))*epsR(startInd-1);
            else
                epsR(startInd:endInd) = interp1([startInd-1,endInd+1],[epsR(startInd-1),epsR(endInd+1)],startInd:endInd);
            end
        end
    end
end

