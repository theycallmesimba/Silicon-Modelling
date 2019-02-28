function analysis = analyzeEffectiveShuttlingResults(sparams, rhos, Hams)
%ANALYZEEFFECTIVESHUTTLINGRESULTS Summary of this function goes here
%   Detailed explanation goes here
    analysis.rhos = rhos;
    
    [~,~,nDataPoints] = size(rhos);
    
    analysis.fidelity = zeros(1,sparams.nStoreDataFrames);
    
    analysis.orbPur = zeros(1,sparams.nStoreDataFrames);
    analysis.valPur = zeros(1,sparams.nStoreDataFrames);
    analysis.spinPur = zeros(1,sparams.nStoreDataFrames);
    analysis.totPur = zeros(1,sparams.nStoreDataFrames);
    
    analysis.orbExp = zeros(3,sparams.nStoreDataFrames);
    analysis.valExp = zeros(3,sparams.nStoreDataFrames);
    analysis.spinExp = zeros(3,sparams.nStoreDataFrames);
    
    for ii = 1:nDataPoints
        [analysis.orbPur(:,ii), analysis.valPur(:,ii), analysis.spinPur(:,ii), ...
            analysis.totPur(:,ii)] = getPurityValues(sparams, squeeze(rhos(:,:,ii)));
        
        [analysis.orbExp(:,ii), analysis.valExp(:,ii), analysis.spinExp(:,ii)] =...
            getExpectationValues(sparams, squeeze(rhos(:,:,ii)));
        
        analysis.fidelity(:,ii) = getFidelity(sparams, squeeze(rhos(:,:,ii)), squeeze(Hams(:,:,ii)));
    end
end

