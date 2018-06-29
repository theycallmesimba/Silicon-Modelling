function gPulse = getInterpolatedPulseValues( sparams, tVec, interpolants )
%GETBATCHOFINTERPPOTSANDPULSE Summary of this function goes here
%   Detailed explanation goes here

    if strcmp(sparams.interpType,'linear')
        % Get the voltage gate values for the current time index
        gPulse = zeros(sparams.numOfGates,length(tVec));
        for vv = 1:sparams.numOfGates
            gPulse(vv,:) = interpolants{vv}({tVec});
        end
    end
end

