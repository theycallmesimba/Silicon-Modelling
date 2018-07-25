function gPulse = getInterpolatedPulseValues( sparams, qVec, interpolants )
%GETBATCHOFINTERPPOTSANDPULSE Summary of this function goes here
%   Detailed explanation goes here

    if strcmp(sparams.interpType,'linear')
        % Get the voltage gate values for the current time index
        gPulse = zeros(sparams.numOfGates,length(qVec));
        for vv = 1:sparams.numOfGates
            currInterpolant = interpolants{vv};
            gPulse(vv,:) = currInterpolant({qVec});
        end
    end
end

