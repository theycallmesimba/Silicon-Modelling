function [sparams, detuningPulse, optPulseTime] = getDetuningPulseAdiabatic(...
    sparams, adiabThreshVec, detBounds, findAdiabatic, effHamiltonianParams )
%GETDETUNINGPULSEADIABATIC Cousin function to getVoltagePulseAdiabatic but
%instead for finding a pulse based solely on detuning.  To be used when we
%don't want to define a detuning pulse in terms of voltage inputs.  This
%code can only handle shuttling between 2 dots and not arbitrary dot chains
%like the cousin function.
    changedSecondSpinFlag = 0;
    if sparams.includeSecondSpin
        sparams.includeSecondSpin = 0;
        changedSecondSpinFlag = 1;
    end

    % Build the outline desired detuning pulse
    ptsForPulse = sparams.nPulsePoints;
    % Left detuning pulse
    tempL = linspace(detBounds(1,1),detBounds(1,1),ptsForPulse/2);
    tempL = [tempL(1:end-1), linspace(detBounds(1,1),detBounds(1,2),ptsForPulse/2)];
    % Right detuning pulse
    tempR = linspace(detBounds(2,2),detBounds(2,1),ptsForPulse/2);
    tempR = [tempR(1:end-1), linspace(detBounds(2,1),detBounds(2,1),ptsForPulse/2)];
    
    detPulse = [tempL; tempR];
    
    [rows,cols] = size(detPulse);
    qXPoints = linspace(0,1,cols);
    adiabQPoints = ones(1,cols)*adiabThreshVec(1);
    
    fig = figure;
    hold on;
    xlabel('Time [arb]','Interpreter','Latex');
    yyaxis left
    for ii = 1:rows
        plot(qXPoints,detPulse(ii,:)/sparams.ee,'--o','Linewidth',1.5);
    end
    ylabel('Detuning [eV]','Interpreter','Latex');
    yyaxis right
    semilogy(qXPoints,adiabQPoints);
    ylabel('Adiabatic parameter [arb]','Interpreter','Latex');
    drawnow;
    pause(1);
    delete(fig);
        
    % This flag is for if we just want a linear collection of detuning 
    % points and don't want to find the corresponding adiabatic pulse
    if ~findAdiabatic
        detuningPulse = detPulse;
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
    optimalPulseDetuning = zeros(rows,length(qXPoints)+1); % +1 for the end points of the pulse
    optimalPulseDetuning(optimalPulseDetuning == 0) = NaN;
    optimalPulseTime = zeros(1,length(qXPoints)+1); % +1 for the end points of the pulse
    optimalPulseTime(optimalPulseTime == 0) = NaN;
    % Initialize our pulse by setting to the same initial conditions
    optimalPulseDetuning(:,1) = detPulse(:,1);
    optimalPulseTime(1) = 0;
    
    fig = figure;
    hold on;
    set(gca,'Fontsize',14,'TickLabelInterpreter','latex');
    ylabel('Detuning [eV]','Interpreter','Latex','Fontsize',20);
    xlabel('Time [s]','Interpreter','Latex','Fontsize',20);
    animatedLines = gobjects(1,sparams.numOfGates);
    color = {'r','g','b','y','m'};
    for ii = 1:2
        animatedLines(ii) = animatedline(optimalPulseTime(1),optimalPulseDetuning(ii,1)/sparams.ee,...
            'Linestyle','--','Marker','o','Linewidth',1.5,'Color',color{ii});
    end
    drawnow;

    for ii = 1:length(qXPoints)        
        % For each point in the voltage pulse, find what time gives it the
        % adiabatic threshold parameter value
        [optimalTimePower, ~] = optimizeTimeForAdiabicity(sparams,ii,...
            length(qXPoints),detPulse,adiabQPoints,effHamiltonianParams);
        optimalTime(ii) = 10^optimalTimePower;

        % Get the corresponding dV/dt for whatever optimal time we found
        % above at the current index point on the voltage pulse
        tVec = linspace(0,optimalTime(ii),length(qXPoints));
        dt = tVec(2) - tVec(1);
        
        % Need to set a special condition for the end of the pulse
        if ii == length(qXPoints)
            optimalPulseDetuning(:,end) = detPulse(:,end);
            optimalPulseTime(end) = optimalPulseTime(end-1) + dt/2;
            break;
        end
        
        % Now find the midway voltage pulse between the current gPulse
        % point and the desired point
        dDet = detPulse(:,ii+1) - detPulse(:,ii);
        optimalPulseDetuning(:,ii+1) = detPulse(:,ii) + dDet./2;
        optimalPulseTime(ii+1) = optimalPulseTime(ii) + dt/2;
        
        for jj = 1:2
            addpoints(animatedLines(jj),optimalPulseTime(ii+1),optimalPulseDetuning(jj,ii+1)/sparams.ee);
        end
        drawnow;
    end
    delete(fig);
    toc;
    
    detuningPulse = zeros(2,sparams.nPulsePoints);
    for ii = 1:2
        detuningPulse(ii,:) = interp1(optimalPulseTime,...
            optimalPulseDetuning(ii,:),linspace(0,max(optimalPulseTime),sparams.nPulsePoints));
    end
    
    optPulseTime = max(optimalPulseTime);
    
    if changedSecondSpinFlag
        sparams.includeSecondSpin = 1;
    end
end

function [time, adiabParams] = optimizeTimeForAdiabicity(sparams, timeInd, numTimeInd,...
    detPulse, thresholdVector, effHamiltonianParams)
    
    % Get order of magnitude of threshold to set tolerance level
    options = optimset('TolX',10^floor(log10(thresholdVector(timeInd))));

    if strcmp(sparams.adiabaticPulseType,'effective')
        [time, ~] = fminbnd(@(qTimePower) optimizeAdiabaticParameterEffective(qTimePower, sparams.nnIndices),...
            sparams.timePowerBounds(1), sparams.timePowerBounds(2), options);
    end

    function diffWithThreshold = optimizeAdiabaticParameterEffective(qTimePower, nn)
        % First let's make the time vector and get the current time we are
        % interested in
        tVec = linspace(0,10^qTimePower,numTimeInd);
        cTime = tVec(timeInd);

        % Make a set of detuning interpolant objects based on the current
        % time vector
        detInterps = makeDetuningInterpolants(tVec, detPulse(1,:), detPulse(2,:));
        
        [adiabaticParam, adiabParams] = calculateAdiabaticParameterEffective(...
            sparams, detInterps, cTime, nn, effHamiltonianParams );

        diffWithThreshold = abs(adiabaticParam - thresholdVector(timeInd));
    end
end
