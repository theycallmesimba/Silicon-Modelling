function analyzePulseTunnelCoupling( sparams, xx, vPulse, pulseTime )
%ANALYZEPULSETUNNELCOUPLING Summary of this function goes here
%   Detailed explanation goes here

    pulseTVec = linspace(0,pulseTime,length(vPulse(1,:)));
    pulseInterps = makePulseInterpolants(sparams,pulseTVec,vPulse);

    cTimes = linspace(0,pulseTime,125);
    currPulse = getInterpolatedPulseValues(sparams,cTimes,pulseInterps);
    
    
    tcs = zeros(1,length(cTimes));
    eL = zeros(1,length(cTimes));
    eR = zeros(1,length(cTimes));
    for ii = 1:length(cTimes)
        currPot = sparams.P2DEGInterpolant(getInterpolantArgument(currPulse(:,ii),xx));
        currPot = squeezeFast(sparams.numOfGates,currPot)';
        [tcs(ii), eL(ii), eR(ii)] = calculateTunnelCoupling(sparams, xx, currPot);
    end
    tcs = tcs./sparams.ee;
    eL(eL == 0) = NaN;
    eR(eR == 0) = NaN;
    
    figure;
    hold on;
    plot(cTimes,eL./sparams.ee);
    plot(cTimes,eR./sparams.ee);
    
    figure;
    plot(cTimes,abs(eL - eR)./sparams.ee);
    
    figure;
    axisFontSize = 20;
    tickFontSize = 15;
    set(gca,'Fontsize',tickFontSize);
    hold on;
    
    yyaxis left
    plot(cTimes,tcs,'Linewidth',2.0,'Linestyle','-');
    
    yyaxis right
    for ii = 1:sparams.numOfGates
        plt = plot(pulseTVec,vPulse(ii,:),'Linewidth',2.0,...
            'Linestyle','-','Marker','none');
        plt.Color(4) = 0.7;
    end
    ylabel('Voltage [V]','Interpreter','Latex','Fontsize',axisFontSize);
end

