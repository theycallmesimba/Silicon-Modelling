function [sweepVec, K, K2] = initializeShuttlingSimulation( sparams, pp )
%INITIALIZESHUTTLINGSIMULATION Summary of this function goes here
%   Detailed explanation goes here

    % Get vector to sweep over during simulation
    sweepVec = getSweepVector(sparams);
    
    % Get momentum operator for simulation
    K = exp(-1i*sparams.dt/2*(pp.^2)/(2*sparams.me*sparams.hbar));
    K2 = K.*K;
    
    % Make the fidelity array
    sparams.fidelity = zeros(length(sweepVec),sparams.nFidelityFrames);
    
    % Mark Stark shift data arrays
    if sparams.calculateStarkShift
        sparams.starkShift = zeros(length(sweepVec),sparams.nStarkShiftFrames);
        sparams.avgEzGround = zeros(length(sweepVec),sparams.nStarkShiftFrames);
        sparams.avgEz = zeros(length(sweepVec),sparams.nStarkShiftFrames);
        sparams.vShiftGround = zeros(length(sweepVec),sparams.nStarkShiftFrames);
        sparams.vShift = zeros(length(sweepVec),sparams.nStarkShiftFrames);
    end
    
    % Let's create the folder to save data
    time = clock;
    sparams.saveFolder = sprintf('%d-%02d-%02d-%02d-%02d-%02d',time(1),time(2),time(3),time(4),time(5),round(time(6)));
    mkdir([sparams.saveDir sparams.saveFolder]);
    sparams.saveFolder = [sparams.saveFolder '/'];
    
    % If we are not sweeping over time for our shuttling simulations, then
    % we need a vector to store the times for our pulses
    if strcmp(sparams.sweptParameter,'adiabicity')
        sparams.totalTime = zeros(1,length(sweepVec));
    end
end

