function analyzePulseAdiabicity( sparams, xx, vPulse, pulseTime, nn, mmVec )
%ANALYZEPULSEADIABICITY Summary of this function goes here
%   Detailed explanation goes here

    pulseTVec = linspace(0,pulseTime,length(vPulse(1,:)));
    pulseInterps = makePulseInterpolants(sparams,pulseTVec,vPulse);
    
    pulseTVec = linspace(0,pulseTime,1000);
    vPulse = getInterpolatedPulseValues(sparams,...
        pulseTVec,pulseInterps);
    
    
    adiabaticParams = zeros(1,length(pulseTVec));
    
    h = waitbar(0,sprintf('Current Time Index: %d/%d',0,length(pulseTVec)),...
        'Name','Analyzing Pulse Adiabicity...',...
        'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
    
    for ii = 1:length(pulseTVec)
        % Check for cancel button click
        if getappdata(h,'canceling')
            flag = 1;
            break;
        end
        
        waitbar(ii/length(pulseTVec), h, sprintf('Current Time Index: %d/%d',ii,length(pulseTVec)));
        
        [adiabaticParams(ii), ~] = calculateAdiabaticParameter( sparams, xx,... 
            pulseInterps, pulseTVec(ii), nn, mmVec );
    end
    
    % Delete the waitbar
    delete(h);
    
    plotFunctionOverVoltagePulse(sparams, pulseTVec, vPulse, adiabaticParams);
end

