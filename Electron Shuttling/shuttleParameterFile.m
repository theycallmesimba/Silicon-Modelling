% % Constants
sparams.ee = 1.602E-19; % Charge of electron
sparams.me = 9.11E-31*0.191; % Effective mass of electron
sparams.hbar = 6.626E-34/(2*pi); % Reduced Planck's constant
sparams.c = 2.998E8; % Speed of light

% % Miscellaneous simulation parameters
sparams.verbose = 1;
sparams.saveData = 1;
sparams.updateWaitbar = 10000; % How many frames to wait before updating waitbar
sparams.updateFigure = 7000; % How many frames to wait before updating figure
sparams.updateInterpPot = 40000; % Interpolate the potential at 40000 points at a time
sparams.nFidelityFrames = 1000; % How many frames to calculate fidelity
sparams.nFigureFrames = 100; % How many frames of the simulation to save for outputted gif

% % Stark Shift Parameters
sparams.calculateStarkShift = 1; % Set to 0 if you don't want to calcualte Stark shift, 1 if you do.
sparams.nStarkShiftFrames = 500; % How many frames of the simulation to calculate the Stark shift
sparams.v0 = 40E9; % [Hz] Reference frequency for Stark shift calculation
sparams.n2 = 2.2*1E-9*1E-9; % [m^2/V^2] Stark shift formula quadratic coefficient

% % Pulse Parameters
sparams.tcThreshold = 0.0009; % Threshold value for the tunnel coupling
% sparams.voltagesToLoad = [0.6,0.693,0.785,0.8]; % Declare which simulated voltages you wish to load
sparams.voltagesToLoad = 0.7:0.01:0.81; % Declare which simulated voltages you wish to load

% % Adiabatic Parameters
sparams.findAdiabaticPulse = 1; % 1 if you want to use adiabatic conditions to find a pulse
sparams.adiabaticThreshold = [0.01,0.0075,0.005,0.0025,0.001];
sparams.hdt = 1E-15; % The h value used to find the time derivative of the wave function
sparams.timePowerBounds = [-14,-6]; % Range of 10^N time values to search in when doing optimizing the adiabatic parameter

% % Shuttling Parameters
sparams.dt = 5E-17; % Time between each simulation frame [sec]
% Vector of times to do the shuttling simulation over
% sparams.totalTime = 8E-9;

% % Directory Parameters
% Where to load potentials from for simulation
% sparams.potDir = 'C:/Users/bbuonaco/Documents/MATLAB/Simulated Potentials/5-gate sim new/';
% sparams.potDir = 'C:/Users/bbuonaco/Documents/MATLAB/Simulated Potentials/Five Gate Gap Induced Quantum Dots - DSize40_GSize20/';
sparams.potDir = 'C:/Users/bbuonaco/Documents/MATLAB/Simulated Potentials/Three Gate Gap Induced Quantum Dots - DSize40_GSize20/';
sparams.interpPotentials = 1; % Whether or not to interpolate the potentials in time domain for simulation
sparams.interpType = 'linear'; % What type of interpolation of the potentials to do in time domain
sparams.numOfGates = 3;
sparams.gateLabels = {'V_1','V_2','V_3'};%,'V_4','V_5'};
sparams.saveDir = 'C:/Users/bbuonaco/Documents/MATLAB/SimulationResults/Electron Shuttling/';