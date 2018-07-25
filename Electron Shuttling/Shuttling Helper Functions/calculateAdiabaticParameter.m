function adiabaticParam = calculateAdiabaticParameter( sparams, xx,... 
    pulseInterpolants, cTime, nn, mmVec )
    
    h = sparams.hdt;

    % Get current pulse point plus the shifted ones as well (+h,+2h,-h,-2h)
    currPulse = getInterpolatedPulseValues(sparams,...
        [cTime-(2*h),cTime-h,cTime,cTime+h,cTime+(2*h)],pulseInterpolants);
    currPot = sparams.P2DEGInterpolant(getInterpolantArgument(currPulse(:,3),xx));
    currPot = squeezeFast(sparams.numOfGates,currPot)';

    [mWFs,mEns] = solve1DSingleElectronSE(sparams, max(mmVec), xx, currPot);

    % Estimate the derivative of the nth WF using a five point stencil
    currPotp2h = sparams.P2DEGInterpolant(getInterpolantArgument(currPulse(:,5),xx));
    currPotp2h = squeezeFast(sparams.numOfGates,currPotp2h)';
    [fxp2h,~] = solve1DSingleElectronSE(sparams, 1, xx, currPotp2h);

    currPotph = sparams.P2DEGInterpolant(getInterpolantArgument(currPulse(:,4),xx));
    currPotph = squeezeFast(sparams.numOfGates,currPotph)';
    [fxph,~] = solve1DSingleElectronSE(sparams, 1, xx, currPotph);

    currPotmh = sparams.P2DEGInterpolant(getInterpolantArgument(currPulse(:,2),xx));
    currPotmh = squeezeFast(sparams.numOfGates,currPotmh)';
    [fxmh,~] = solve1DSingleElectronSE(sparams, 1, xx, currPotmh);

    currPot2mh = sparams.P2DEGInterpolant(getInterpolantArgument(currPulse(:,1),xx));
    currPot2mh = squeezeFast(sparams.numOfGates,currPot2mh)';
    [fxm2h,~] = solve1DSingleElectronSE(sparams, 1, xx, currPot2mh);

    nWFDeriv = (-fxp2h + 8*fxph - 8*fxmh + fxm2h)/(12*h);

    % Now evaluate the adiabatic parameter by summing over all excited
    % states
    adiabaticParam = 0;
    for jj = mmVec
        if jj == nn
            continue
        end
        adiabaticParamTemp = abs(sparams.hbar*getInnerProduct(xx,mWFs(:,jj),nWFDeriv)/(mEns(jj,jj) - mEns(nn,nn)));
        adiabaticParam = adiabaticParam + adiabaticParamTemp;
    end
end