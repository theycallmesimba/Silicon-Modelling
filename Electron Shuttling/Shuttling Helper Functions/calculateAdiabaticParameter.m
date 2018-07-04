function calculateAdiabaticParameter( sparams, xx )
%FINDADIABATICTHRESHOLD Summary of this function goes here
%   Detailed explanation goes here
    profile on
    
    time = 8E-9;
    h = 1E-14;
    qTime = linspace(0,time,501);
    thresh = 0;
    numOfMs = 5;
    nn = 1;
    
    % Now we need to make the individual gate interpolants for the pulse
    tPots = linspace(0,time,length(sparams.voltagePulse(1,:)));
    vPulseGInterpolants = {};
    for vv = 1:sparams.numOfGates
        vPulseGInterpolants{vv} = griddedInterpolant({tPots},sparams.voltagePulse(vv,:));
    end
    
    adiabaticThresholdValue = zeros(1,length(qTime));
    for ii = 1:length(qTime)
        currPulse = getInterpolatedPulseValues(sparams,qTime(ii),vPulseGInterpolants);
        currPot = sparams.P2DEGInterpolant(getInterpolantArgument(currPulse,xx));
        currPot = squeezeFast(sparams.numOfGates,currPot)';

        [mWFs,mEns] = solve1DSingleElectronSE(sparams, numOfMs, xx, currPot);

        % Estimate the derivative of the nth WF using a five point stencil
        currPulse = getInterpolatedPulseValues(sparams,qTime(ii)+(2*h),vPulseGInterpolants);
        currPotH = sparams.P2DEGInterpolant(getInterpolantArgument(currPulse,xx));
        currPotH = squeezeFast(sparams.numOfGates,currPotH)';
        [fxp2h,~] = solve1DSingleElectronSE(sparams, 1, xx, currPotH);

        currPulse = getInterpolatedPulseValues(sparams,qTime(ii)+(h),vPulseGInterpolants);
        currPotH = sparams.P2DEGInterpolant(getInterpolantArgument(currPulse,xx));
        currPotH = squeezeFast(sparams.numOfGates,currPotH)';
        [fxph,~] = solve1DSingleElectronSE(sparams, 1, xx, currPotH);

        currPulse = getInterpolatedPulseValues(sparams,qTime(ii)-(h),vPulseGInterpolants);
        currPotH = sparams.P2DEGInterpolant(getInterpolantArgument(currPulse,xx));
        currPotH = squeezeFast(sparams.numOfGates,currPotH)';
        [fxmh,~] = solve1DSingleElectronSE(sparams, 1, xx, currPotH);

        currPulse = getInterpolatedPulseValues(sparams,qTime(ii)-(2*h),vPulseGInterpolants);
        currPotH = sparams.P2DEGInterpolant(getInterpolantArgument(currPulse,xx));
        currPotH = squeezeFast(sparams.numOfGates,currPotH)';
        [fxm2h,~] = solve1DSingleElectronSE(sparams, 1, xx, currPotH);

        nWFDeriv = (-fxp2h + 8*fxph - 8*fxmh + fxm2h)/(12*h);

        thresh = 0;
        for jj = 1:numOfMs
            if jj == nn
                continue
            end
            threshTemp = abs(sparams.hbar*getInnerProduct(xx,mWFs(:,jj),nWFDeriv)/(mEns(jj,jj) - mEns(nn,nn)));
            thresh = thresh + threshTemp;
        end
        adiabaticThresholdValue(ii) = thresh;
    end
    
    plotFunctionOverVoltagePulse(sparams,tPots,qTime,adiabaticThresholdValue);
    
    profile off
    profile viewer
end






