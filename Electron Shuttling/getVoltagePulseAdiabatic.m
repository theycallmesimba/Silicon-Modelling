function [sparams, optPulseTime] = getVoltagePulseAdiabatic( sparams, xx, adiabaticThreshVec )
%GETVOLTAGEPULSE Summary of this function goes here
%   Detailed explanation goes here

    % This function is to generate a pulsing sequence for
    % all the gates we have control over in our geometry.    
    sparams.adiabaticThresholdValue = 0.001;    
    sparams.tcMax = zeros(sparams.numOfGates-1,1);
    
    gMin = ones(1,sparams.numOfGates)*0.6;
    gMax = zeros(1,sparams.numOfGates);
    
    gMax(1) = 0.800;
    vMinBnd = 0.795;
    vMaxBnd = 0.810;
    
    % Find the gate voltage that maximizes the tunnel coupling for each
    % gate
    for ii = 2:sparams.numOfGates
        gArg = gMin;
        gArg(ii-1) = gMax(ii-1);
        [gMax(ii), sparams.tcMax(ii-1)] = findTunnelCouplingMax(sparams, xx, gArg, vMinBnd, vMaxBnd, ii);
    end
    
    % Now we will construct the rough outline of the pulse with the main
    % points that we require the voltage pulse to sweep through
    [gPulse, xPoints] = formRoughPulseOutline(sparams, gMin, gMax);
    
    % Now we want to interpolate our pulse a bit more to have more points
    % to calculate the adiabatic parameter at
    pulseInterps = makePulseInterpolants(sparams,xPoints,gPulse);
    fracBuf = 0.035;
    fracDens = 0.35;
    ptsPerGate = 100;
    
    ptsBuf = xPoints(2)*fracBuf;
    nPtsSparse = round(ptsPerGate/2*fracDens);
    nPtsDense = round(ptsPerGate/2) - nPtsSparse;
    
    qXPoints = linspace(0,xPoints(2)-ptsBuf,nPtsSparse);
%     whichThreshold = ones(1,nPtsSparse)*adiabaticThreshVec(1);
    qXPoints = [qXPoints, linspace(xPoints(2)-ptsBuf,...
        xPoints(2),nPtsDense)];
    
    whichThreshold = logspace(log10(adiabaticThreshVec(1)),log10(adiabaticThreshVec(2)),nPtsDense + nPtsSparse);
%     whichThreshold = [whichThreshold, ones(1,nPtsDense)*adiabaticThreshVec(1)];
    for ii = 1:(length(xPoints)-3)/2
        qXPoints = [qXPoints, linspace(xPoints(2*ii),...
            xPoints(2*ii)+ptsBuf,nPtsDense)];
%         whichThreshold = [whichThreshold, ones(1,nPtsDense)*adiabaticThreshVec(1)];
        
        qXPoints = [qXPoints, linspace(xPoints(2*ii)+ptsBuf,...
            xPoints(2*ii+1),nPtsSparse)];
        whichThreshold = [whichThreshold...
            logspace(log10(adiabaticThreshVec(2)),log10(adiabaticThreshVec(1)),nPtsDense + nPtsSparse)];
        
        qXPoints = [qXPoints, linspace(xPoints(2*ii+1),...
            xPoints(2*(ii+1))-ptsBuf,nPtsSparse)];
%         whichThreshold = [whichThreshold, ones(1,2*nPtsSparse)*adiabaticThreshVec(1)];
        qXPoints = [qXPoints, linspace(xPoints(2*(ii+1))-ptsBuf,...
            xPoints(2*(ii+1)),nPtsDense)];
        whichThreshold = [whichThreshold...
            logspace(log10(adiabaticThreshVec(1)),log10(adiabaticThreshVec(2)),nPtsDense + nPtsSparse)];
%         whichThreshold = [whichThreshold, ones(1,nPtsDense)*adiabaticThreshVec(1)];
    end
    qXPoints = [qXPoints, linspace(xPoints(end-1),...
            xPoints(end-1)+ptsBuf,nPtsDense)];
%     whichThreshold = [whichThreshold, ones(1,nPtsDense)*adiabaticThreshVec(1)];
    qXPoints = [qXPoints, linspace(xPoints(end-1)+ptsBuf,...
            xPoints(end),nPtsSparse)];
    whichThreshold = [whichThreshold...
            logspace(log10(adiabaticThreshVec(2)),log10(adiabaticThreshVec(1)),nPtsDense + nPtsSparse)];
%     whichThreshold = [whichThreshold, ones(1,nPtsSparse)*adiabaticThreshVec(1)];
        
    gPulse = getInterpolatedPulseValues(sparams, qXPoints, pulseInterps);
    
    fig = figure;
    hold on;
    yyaxis left
    for ii = 1:sparams.numOfGates
        plot(qXPoints,gPulse(ii,:),'--o','Linewidth',1.5);
    end
    yyaxis right
    plot(qXPoints,whichThreshold);
    drawnow;
    pause(10);
    delete(fig);
    
    % Now, we want to go through each point in our pulse and find what the
    % derivatve of each individual gate pulse should be at each point in
    % time.  We do this by fixing an h (i.e. setting our dt to be fixed)
    % and vary what f(t) looks like by adjusting the dV for each gate
    % voltage thereby obtaining dV/dt at each point in the pulse.
    tic;
    optimalTime = zeros(1,length(qXPoints));
    optimalPulseVoltage = zeros(sparams.numOfGates,length(qXPoints)+1); % +1 for the end points of the pulse
    optimalPulseVoltage(optimalPulseVoltage == 0) = NaN;
    optimalPulseTime = zeros(1,length(qXPoints)+1); % +1 for the end points of the pulse
    optimalPulseTime(optimalPulseTime == 0) = NaN;
    % Initialize our pulse by setting to the same initial conditions
    optimalPulseVoltage(:,1) = gPulse(:,1);
    optimalPulseTime(1) = 0;
    
    fig = figure;
    hold on;
    ylabel('Voltage [V]','Interpreter','Latex');
    xlabel('Time [s]','Interpreter','Latex');
    animatedLines = gobjects(1,sparams.numOfGates);
    color = {'r','g','b','y','m'};
    for ii = 1:sparams.numOfGates
        animatedLines(ii) = animatedline(optimalPulseTime(1),optimalPulseVoltage(ii,1),...
            'Linestyle','--','Marker','o','Linewidth',1.5,'Color',color{ii});
    end
    drawnow;
    
    for ii = 1:length(qXPoints)        
        % For each point in the voltage pulse, find what time gives it the
        % adiabatic threshold parameter value
        optimalTimePower = optimizeTimeForAdiabicity(sparams,xx,ii,...
            length(qXPoints),gPulse,whichThreshold);
        optimalTime(ii) = 10^optimalTimePower;

        % Get the corresponding dV/dt for whatever optimal time we found
        % above at the current index point on the voltage pulse
        tVec = linspace(0,optimalTime(ii),length(qXPoints));
        dt = tVec(2) - tVec(1);
        dV = zeros(sparams.numOfGates,1);
        
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
        
        for jj = 1:sparams.numOfGates
            addpoints(animatedLines(jj),optimalPulseTime(ii+1),optimalPulseVoltage(jj,ii+1));
        end
        drawnow;
    end
    delete(fig);
    toc;
    
    currPulsePtsNum = length(optimalPulseVoltage(1,:));
    sparams.voltagePulse = zeros(sparams.numOfGates,2*currPulsePtsNum);
    for ii = 1:sparams.numOfGates
        sparams.voltagePulse(ii,:) = interp1(optimalPulseTime,...
            optimalPulseVoltage(ii,:),linspace(0,max(optimalPulseTime),length(sparams.voltagePulse(1,:))));
    end
    
    optPulseTime = max(optimalPulseTime);
end

function time = optimizeTimeForAdiabicity(sparams, xx, timeInd, numTimeInd,...
    gPulse, thresholdVector)
    
    % Get order of magnitude of threshold to set tolerance level
    options = optimset('TolX',10^floor(log10(thresholdVector(timeInd))));
    
    [time, ~] = fminbnd(@(qTimePower) optimizeAdiabaticParameter(qTimePower, 1, [2,3,4]),...
        sparams.timePowerBounds(1), sparams.timePowerBounds(2), options);

    function diffWithThreshold = optimizeAdiabaticParameter(qTimePower, nn, mmVec)
        % First let's make the time vector and get the current time we are
        % interested in
        tVec = linspace(0,10^qTimePower,numTimeInd);
        cTime = tVec(timeInd);
        
        % Make a new set of pulse interpolant objects based on the current time
        % vector
        pulseInterps = makePulseInterpolants(sparams,tVec,gPulse);
        
        adiabaticParam = calculateAdiabaticParameter( sparams, xx,...
            pulseInterps, cTime, nn, mmVec );
        
        diffWithThreshold = abs(adiabaticParam - thresholdVector(timeInd));
    end
end

function [gPulse, xPoints] = formRoughPulseOutline(sparams, gMin, gMax)
    gPulse = [];
    for ii = 1:sparams.numOfGates
        % Pulse point where only one gate voltage is high
        gArg = gMin;
        gArg(ii) = gMax(ii);
        gPulse = [gPulse, gArg'];
        
        if ii == sparams.numOfGates
            break;
        end
        
        % Pulse point where two gate voltages are high and tc is maxed
        gArg(ii+1) = gMax(ii+1);
        gPulse = [gPulse, gArg'];
    end
    
    [~,cols] = size(gPulse);
    xPoints = linspace(0,1,cols);
    
    fig = figure;
    hold on;
    for ii = 1:sparams.numOfGates
        plot(xPoints,gPulse(ii,:),'--o','Linewidth',1.5);
    end
    drawnow;
    pause(2.5);
    delete(fig);
end
