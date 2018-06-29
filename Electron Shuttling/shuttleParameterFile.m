% % Constants
sparams.ee = 1.602E-19; % Charge of electron
sparams.me = 9.11E-31*0.191; % Effective mass of electron
sparams.hbar = 6.626E-34/(2*pi); % Reduced Planck's constant
sparams.c = 2.998E8; % Speed of light

sparams.v0 = 40E9; % [Hz] Reference frequency for Stark shift calculation
sparams.n2 = 2.2*1E-9*1E-9; % [m^2/V^2] Stark shift formula quadratic coefficient taken from Dzurak's paper

sparams.verbose = 1;

% % Miscellaneous simulation parameters
sparams.updateWaitbar = 10000; % How many frames to wait before updating waitbar
sparams.updateFigure = 7000; % How many frames to wait before updating figure
sparams.updateInterpPot = 40000; % Interpolate the potential at 40000 points at a time
sparams.nFidelityFrames = 1000; % How many frames to calculate fidelity
sparams.nFigureFrames = 160; % How many frames of the simulation to save for outputted gif
sparams.nStarkShiftFrames = 500; % How many frames of the simulation to calculate the Stark shift
sparams.calculateStarkShift = 1; % Set to 0 if you don't want to calcualte Stark shift, 1 if you do.

sparams.tcTuningThreshold = 0.0011;

% % Main simulation parameters
% sparams.potFile = 'Shuttling3gates_121steps.xlsx'; % Files to load potentials from for simulation
% sparams.potDir = 'simulatedPotentials/0427/'; % Directory where the pot files are
% sparams.potDir = 'C:/Users/bbuonaco/Documents/MATLAB/Simulated Potentials/Five Gate Gap Induced Quantum Dots - DSize40_GSize20/';
sparams.potDir = 'C:/Users/bbuonaco/Documents/MATLAB/Simulated Potentials/Five Gate Gap Induced Quantum Dots - DSize40_GSize20/';
sparams.interpPotentials = 1; % Whether or not to interpolate the potentials in time domain for simulation
sparams.interpType = 'linear'; % What type of interpolation of the potentials to do in time domain
sparams.numOfGates = 5;
sparams.gateLabels = {'V_{1}','V_{2}','V_{3}','V_{4}','V_{5}'};

sparams.saveDir = 'C:/Users/bbuonaco/Documents/MATLAB/SimulationResults/Electron Shuttling/';
sparams.saveData = 1;

sparams.dt = 5E-17; % Time between each simulation frame [sec]
% sparams.totalTime = logspace(-11,-8.15,25); % Total time of simulted pulse [sec]
% sparams.totalTime = [1E-11,2E-11];
%sparams.totalTime = 3E-9;
sparams.totalTime = 8E-9;
% sparams.totalTime = [4.5973E-10,4.0973E-9];