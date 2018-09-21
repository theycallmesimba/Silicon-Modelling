function sweepVec = getSweepVector( sparams )
%GETSWEEPVECTOR Summary of this function goes here
%   Detailed explanation goes here

    if strcmp(sparams.sweptParameter,'time') 
        sweepVec = sparams.totalTime;
    elseif strcmp(sparams.sweptParameter,'adiabicity')
        sweepVec = sparams.adiabaticThreshold;
    end
end

