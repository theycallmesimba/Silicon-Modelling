function [sparams, voltagePulse, optPulseTime] = getVoltagePulseAdiabatic(...
    sparams, xx, adiabThreshVec, vBounds, findAdiabatic, effHamiltonianParams )
%GETVOLTAGEPULSE Summary of this function goes here
%   Detailed explanation goes here
    
    try
        changedSecondSpinFlag = 0;
        if sparams.includeSecondSpin
            sparams.includeSecondSpin = 0;
            changedSecondSpinFlag = 1;
        end
    catch
        % Just ignore it because it means we aren't doing an effective
        % shuttling simulation
    end
    
    % This function is to generate a pulsing sequence for
    % all the gates we have control over in our geometry.    
    nGatesUsedInPulse = length(sparams.gatesUsedInPulse);
    sparams.tcMax = zeros(nGatesUsedInPulse-1,1);
    
    gMin = ones(1,nGatesUsedInPulse);
    for ii = 1:nGatesUsedInPulse
        gMin(1,ii) = vBounds(ii,1);
    end
    gMax = gMin;
    
    gMax(sparams.gatesUsedInPulse(1)) = vBounds(1,2);
    
    % Find the gate voltage that maximizes the tunnel coupling for each
    % gate
    for ii = 2:nGatesUsedInPulse
        gArg = gMin;
        gArg(sparams.gatesUsedInPulse(ii-1)) = gMax(sparams.gatesUsedInPulse(ii-1));

        [gMax(sparams.gatesUsedInPulse(ii)), sparams.tcMax(ii-1)] = findResonantTunnelCoupling(sparams, xx, gArg, sparams.gatesUsedInPulse(ii), 1);
    end
    
    % Now we will construct the rough outline of the pulse with the main
    % points that we require the voltage pulse to sweep through
    [gPulse, xPoints] = formRoughPulseOutline(sparams, gMin, gMax);
    
    % Now we want to interpolate our pulse a bit more to have more points
    % to calculate the adiabatic parameter at
    pulseInterps = makePulseInterpolants(sparams,xPoints,gPulse);
    fracBuf = 0.035;
    fracNPoints = 0.35;
    ptsPerGate = 150;
    
    ptsBuf = xPoints(2)*fracBuf;
    nPtsSparse = round(ptsPerGate/2*fracNPoints);
    nPtsDense = round(ptsPerGate/2) - nPtsSparse;
    
    qXPoints = linspace(0,xPoints(2)-ptsBuf,nPtsSparse);
    qXPoints = [qXPoints, linspace(xPoints(2)-ptsBuf,xPoints(2),nPtsDense)];
    for ii = 1:(length(xPoints)-3)/2
        qXPoints = [qXPoints, linspace(xPoints(2*ii),...
            xPoints(2*ii)+ptsBuf,nPtsDense)];
        
        qXPoints = [qXPoints, linspace(xPoints(2*ii)+ptsBuf,...
            xPoints(2*ii+1),nPtsSparse)];
        
        qXPoints = [qXPoints, linspace(xPoints(2*ii+1),...
            xPoints(2*(ii+1))-ptsBuf,nPtsSparse)];
        qXPoints = [qXPoints, linspace(xPoints(2*(ii+1))-ptsBuf,...
            xPoints(2*(ii+1)),nPtsDense)];
    end
    qXPoints = [qXPoints, linspace(xPoints(end-1),xPoints(end-1)+ptsBuf,nPtsDense)];
    qXPoints = [qXPoints, linspace(xPoints(end-1)+ptsBuf,xPoints(end),nPtsSparse)];
        
    gPulse = getInterpolatedPulseValues(sparams, qXPoints, pulseInterps);
    
    % Get the vector of points corresponding to what the adiabatic
    % parameter should be at each point in time
    adiabPoints = [];
    for ii = 2:nGatesUsedInPulse
        adiabPoints = [adiabPoints,linspace(adiabThreshVec(1),adiabThreshVec(1),nPtsSparse)];
        adiabPoints = [adiabPoints,logspace(log10(adiabThreshVec(1)),log10(adiabThreshVec(2)),nPtsDense),...
            logspace(log10(adiabThreshVec(2)),log10(adiabThreshVec(1)),nPtsDense)];
        adiabPoints = [adiabPoints,linspace(adiabThreshVec(1),adiabThreshVec(1),nPtsSparse)];
    end
    adiabQPoints = adiabPoints;
    
    fig = figure;
    hold on;
    xlabel('Time [arb]','Interpreter','Latex');
    yyaxis left
    for ii = 1:length(sparams.gatesUsedInPulse)
        plot(qXPoints,gPulse(ii,:),'--o','Linewidth',1.5);
    end
    xlabel('Voltage [V]','Interpreter','Latex');
    yyaxis right
    semilogy(qXPoints,adiabQPoints);
    ylabel('Adiabatic parameter [arb]','Interpreter','Latex');
    drawnow;
    pause(0.75);
    delete(fig);
    
    % This flag is for if we just want a linear collection of voltage 
    % points and don't want to find the corresponding adiabatic pulse
    if ~findAdiabatic
        voltagePulse = gPulse;
        optPulseTime = NaN;
        
        return;
    end
    
    % Now, we want to go through each point in our pulse and find what the
    % derivatve of each individual gate pulse should be at each point in
    % time.  We do this by fixing an h (i.e. setting our dt to be fixed)
    % and vary what f(t) looks like by adjusting the dV for each gate
    % voltage thereby obtaining dV/dt at each point in the pulse.
    tic;
    optimalTime = zeros(1,length(qXPoints));
    optimalPulseVoltage = zeros(length(sparams.gatesUsedInPulse),length(qXPoints)+1); % +1 for the end points of the pulse
    optimalPulseVoltage(optimalPulseVoltage == 0) = NaN;
    optimalPulseTime = zeros(1,length(qXPoints)+1); % +1 for the end points of the pulse
    optimalPulseTime(optimalPulseTime == 0) = NaN;
    % Initialize our pulse by setting to the same initial conditions
    optimalPulseVoltage(:,1) = gPulse(:,1);
    optimalPulseTime(1) = 0;
    
    fig = figure;
    set(gca,'TickLabelInterpreter','latex','Fontsize',14);
    hold on;
    ylabel('Voltage [V]','Interpreter','Latex','Fontsize',22);
    xlabel('Time [s]','Interpreter','Latex','Fontsize',22);
    animatedLines = gobjects(1,length(sparams.gatesUsedInPulse));
    color = {'r','g','b','y','m'};
    for ii = 1:length(sparams.gatesUsedInPulse)
        animatedLines(ii) = animatedline(optimalPulseTime(1),optimalPulseVoltage(ii,1),...
            'Linestyle','--','Marker','o','Linewidth',1.5,'Color',color{ii},'DisplayName',sprintf('$V_%d$',ii));
    end
    leg = legend;
    leg.Interpreter = 'latex';
    drawnow;
    
    for ii = 1:length(qXPoints)        
        % For each point in the voltage pulse, find what time gives it the
        % desired adiabatic threshold parameter value
        [optimalTimePower, ~] = optimizeTimeForAdiabicity(sparams,xx,ii,...
            length(qXPoints),gPulse,adiabPoints,effHamiltonianParams);
        optimalTime(ii) = 10^optimalTimePower;

        % Get the corresponding dV/dt for whatever optimal time we found
        % above at the current index point on the voltage pulse
        tVec = linspace(0,optimalTime(ii),length(qXPoints));
        dt = tVec(2) - tVec(1);
        
        % Need to set a special condition for the end of the pulse
        if ii == length(qXPoints)
            optimalPulseVoltage(:,end) = gPulse(:,end);
            optimalPulseTime(end) = optimalPulseTime(end-1) + dt/2;
            break;
        end
        
        % Now find the midway voltage pulse between the current gPulse
        % point and the desired point
        dV = gPulse(:,ii+1) - gPulse(:,ii);
        optimalPulseVoltage(:,ii+1) = gPulse(:,ii) + dV./2;
        optimalPulseTime(ii+1) = optimalPulseTime(ii) + dt/2;
        
        for jj = 1:length(sparams.gatesUsedInPulse)
            addpoints(animatedLines(jj),optimalPulseTime(ii+1),optimalPulseVoltage(jj,ii+1));
        end
        drawnow;
    end
    delete(fig);
    toc;
    
    voltagePulse = zeros(length(sparams.gatesUsedInPulse),sparams.nPulsePoints);
    for ii = 1:length(sparams.gatesUsedInPulse)
        voltagePulse(ii,:) = interp1(optimalPulseTime,...
            optimalPulseVoltage(ii,:),linspace(0,max(optimalPulseTime),sparams.nPulsePoints));
    end
    
    optPulseTime = max(optimalPulseTime);
    
    if changedSecondSpinFlag
        sparams.includeSecondSpin = 1;
    end
end

function [time, adiabParams] = optimizeTimeForAdiabicity(sparams, xx, timeInd, numTimeInd,...
    gPulse, thresholdVector, effHamiltonianParams)
    
    % Get order of magnitude of threshold to set tolerance level
    options = optimset('TolX',10^floor(log10(thresholdVector(timeInd))));

    if strcmp(sparams.adiabaticPulseType,'coherent') 
        [time, ~] = fminbnd(@(qTimePower) optimizeAdiabaticParameter(qTimePower, 1, 1:4),...
            sparams.timePowerBounds(1), sparams.timePowerBounds(2), options);
    elseif strcmp(sparams.adiabaticPulseType,'effective')
        [time, ~] = fminbnd(@(qTimePower) optimizeAdiabaticParameterEffective(qTimePower, 1),...
            sparams.timePowerBounds(1), sparams.timePowerBounds(2), options);
    end
    
    function diffWithThreshold = optimizeAdiabaticParameter(qTimePower, nn, mmVec)
        % First let's make the time vector and get the current time we are
        % interested in
        tVec = linspace(0,10^qTimePower,numTimeInd);
        cTime = tVec(timeInd);

        % Make a new set of pulse interpolant objects based on the current time
        % vector
        pulseInterps = makePulseInterpolants(sparams,tVec,gPulse);
        
        [adiabaticParam, adiabParams] = calculateAdiabaticParameter(...
            sparams, xx, pulseInterps, cTime, nn, mmVec );

        diffWithThreshold = abs(adiabaticParam - thresholdVector(timeInd));
    end

    function diffWithThreshold = optimizeAdiabaticParameterEffective(qTimePower, nn)
        % First let's make the time vector and get the current time we are
        % interested in
        tVec = linspace(0,10^qTimePower,numTimeInd);
        cTime = tVec(timeInd);

        % Make a set of detuning interpolant objects based on the current
        % time vector
        [epsL, epsR, ~, ~, ~, ~, ~, ~] = decodeEffHamiltonianParamVariable(effHamiltonianParams);
        detInterps = makeDetuningInterpolants(tVec, epsL, epsR);
        
        [adiabaticParam, adiabParams] = calculateAdiabaticParameterEffective(...
            sparams, detInterps, cTime, nn, effHamiltonianParams );

        diffWithThreshold = abs(adiabaticParam - thresholdVector(timeInd));
    end
end

function [gPulse, xPoints] = formRoughPulseOutline(sparams, gMin, gMax)
    gPulse = [];
    
    nGatesUsedInPulse = length(sparams.gatesUsedInPulse);
    for ii = 1:nGatesUsedInPulse
        % Pulse point where only one gate voltage is high
        gArg = gMin;
        gArg(sparams.gatesUsedInPulse(ii)) = gMax(sparams.gatesUsedInPulse(ii));
        gPulse = [gPulse, gArg'];
        
        if ii == nGatesUsedInPulse
            break;
        end
        
        % Pulse point where two gate voltages are high and tc is maxed (i.e
        % 0 detuning)
        gArg(sparams.gatesUsedInPulse(ii+1)) = gMax(sparams.gatesUsedInPulse(ii+1));
        gPulse = [gPulse, gArg'];
    end
    
    [~,cols] = size(gPulse);
    xPoints = linspace(0,1,cols);
    
    fig = figure;
    hold on;
    for ii = 1:nGatesUsedInPulse
        plot(xPoints,gPulse(ii,:),'--o','Linewidth',1.5);
    end
    drawnow;
    pause(0.75);
    delete(fig);
end
