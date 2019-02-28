function gPulse = getInterpolatedPulseValues( sparams, qVec, interpolants )
%GETBATCHOFINTERPPOTSANDPULSE Summary of this function goes here
%   Detailed explanation goes here

    % Get the voltage gate values for the current time index
    gPulse = zeros(length(sparams.gatesUsedInPulse),length(qVec));
    for vv = 1:length(sparams.gatesUsedInPulse)
        currInterpolant = interpolants{vv};
        gPulse(vv,:) = currInterpolant({qVec});
    end
end

