function [adiabParamNNTotal, adiabParamInd] = calculateAdiabaticParameter( sparams, xx,... 
    pulseInterpolants, cTime, nn, mmVec )
    
    h = sparams.hdt;

    % Get current pulse point plus the shifted ones as well (+h,+2h,-h,-2h)
    currPulse = getInterpolatedPulseValues(sparams,...
        [cTime-(2*h),cTime-h,cTime,cTime+h,cTime+(2*h)],pulseInterpolants);
    
    currPot = sparams.P2DEGInterpolant(getInterpolantArgument(currPulse(:,3),xx));
    currPot = squeezeFast(sparams.numOfGates,currPot)';
    [WFs,ens] = solve1DSingleElectronSE(sparams, max(mmVec), xx, currPot);
    
    % Estimate the derivative of the nth WF using a five point stencil
    currPotp2h = sparams.P2DEGInterpolant(getInterpolantArgument(currPulse(:,5),xx));
    currPotp2h = squeezeFast(sparams.numOfGates,currPotp2h)';
    [fxp2h,~] = solve1DSingleElectronSE(sparams, max(mmVec), xx, currPotp2h);

    currPotph = sparams.P2DEGInterpolant(getInterpolantArgument(currPulse(:,4),xx));
    currPotph = squeezeFast(sparams.numOfGates,currPotph)';
    [fxph,~] = solve1DSingleElectronSE(sparams, max(mmVec), xx, currPotph);

    currPotmh = sparams.P2DEGInterpolant(getInterpolantArgument(currPulse(:,2),xx));
    currPotmh = squeezeFast(sparams.numOfGates,currPotmh)';
    [fxmh,~] = solve1DSingleElectronSE(sparams, max(mmVec), xx, currPotmh);

    currPot2mh = sparams.P2DEGInterpolant(getInterpolantArgument(currPulse(:,1),xx));
    currPot2mh = squeezeFast(sparams.numOfGates,currPot2mh)';
    [fxm2h,~] = solve1DSingleElectronSE(sparams, max(mmVec), xx, currPot2mh);
    
    % Now, the SE equation solver uses eigs which can output a random +/-1
    % global phase on the wavefunction output depending on what random seed
    % vector it uses to solve the eigenvalue problem.  We need to make sure
    % our phases are consistent across all dts we calculate otherwise we'll
    % get big jumps in the derivative.  We can check the phase by taking
    % the inner product which (if dt is small) should be ~1.  So if it's <0
    % then we know the phase is off.
    for ii = mmVec
        if getInnerProduct(xx, WFs(:,ii), fxp2h(:,ii)) < 0
            fxp2h(:,ii) = -fxp2h(:,ii);
        end
        if getInnerProduct(xx, WFs(:,ii), fxph(:,ii)) < 0
            fxph(:,ii) = -fxph(:,ii);
        end
        if getInnerProduct(xx, WFs(:,ii), fxmh(:,ii)) < 0
            fxmh(:,ii) = -fxmh(:,ii);
        end
        if getInnerProduct(xx, WFs(:,ii), fxm2h(:,ii)) < 0
            fxm2h(:,ii) = -fxm2h(:,ii);
        end
    end
    
    
    % Now evaluate the adiabatic parameter by summing over all excited
    % states
    adiabParamNNTotal = 0;
    vv = 0;
    adiabParamInd = struct([]);
    % <jj|\dot{ii}> adiabatic parameter
    for ii = mmVec
        % Get |\dot{ii}>
        iiWFDeriv = (-fxp2h(:,ii) + 8*fxph(:,ii) - 8*fxmh(:,ii) + fxm2h(:,ii))/(12*h);
        for jj = mmVec
            % Ignore <ii|\dot{ii}> terms
            if ii == jj
                continue
            end
            vv = vv + 1;
            adiabParamInd(vv).mm = jj;
            adiabParamInd(vv).nn = ii;
            adiabParamInd(vv).adiabaticParameter = sparams.hbar*abs(getInnerProduct(xx,WFs(:,jj),iiWFDeriv)/(ens(ii,ii) - ens(jj,jj)));
            
            if ii == nn
                adiabParamNNTotal = adiabParamNNTotal + adiabParamInd(vv).adiabaticParameter;
            end
        end
    end
end