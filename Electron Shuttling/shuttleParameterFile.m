% % Constants
sparams.ee = 1.602E-19; % Charge of electron
sparams.me = 9.11E-31*0.191; % Effective mass of electron
sparams.hbar = 6.626E-34/(2*pi); % Reduced Planck's constant
sparams.c = 2.998E8; % Speed of light

sparams.v0 = 40E9; % [Hz] Reference frequency for Stark shift calculation
sparams.n2 = 2.2*1E-9*1E-9; % [m^2/V^2] Stark shift formula quadratic coefficient taken from Dzurak's paper

sparams.verbose = 1;

% % Miscellaneous simulation parameters
sparams.updateWaitbar = 5000; % How many frames to wait before updating waitbar
sparams.updateFigure = 1000; % How many frames to wait before updating figure
sparams.updateFidelity = 4000; % How many frames to wait before calculating fidelity
sparams.updateInterpPot = 40000; % Interpolate the potential at 40000 points at a time
sparams.nFigureFrames = 75; % How many frames of the simulation to save for outputted gif
sparams.nStarkShiftFrames = 300; % How many frames of the simulation to calculate the Stark shift
sparams.calculateStarkShift = 1; % Set to 0 if you don't want to calcualte Stark shift, 1 if you do.

% % Main simulation parameters
% sparams.potFile = 'Shuttling3gates_121steps.xlsx'; % Files to load potentials from for simulation
% sparams.potDir = 'simulatedPotentials/0427/'; % Directory where the pot files are
sparams.potDir = '../../Simulated Potentials/Three Gate Gap Induced Quantum Dots - DSize40_GSize20/';
sparams.interpPotentials = 1; % Whether or not to interpolate the potentials in time domain for simulation
sparams.interpType = 'linear'; % What type of interpolation of the potentials to do in time domain
sparams.numOfGates = 3;

sparams.dt = 5E-17; % Time between each simulation frame [sec]
% sparams.totalTime = logspace(-10,-8,10); % Total time of simulted pulse [sec]
% sparams.totalTime = [1E-12,2E-12];
% sparams.totalTime = [1E-9,5E-9,1E-10];
sparams.totalTime = 1E-12;