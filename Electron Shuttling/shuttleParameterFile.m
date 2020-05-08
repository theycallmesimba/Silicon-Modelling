% % Constants
sparams.ee = 1.602E-19; % Charge of electron
sparams.me = 9.11E-31*0.191; % Effective mass of electron
sparams.hbar = 6.626E-34/(2*pi); % Reduced Planck's constant
sparams.c = 2.998E8; % Speed of light
sparams.g = 2; % Electron g-factor
sparams.mub = 9.274E-24; % Bohr's magneton
% Pauli matrices
sparams.sigmaz = [1,0;0,-1];
sparams.sigmax = [0,1;1,0];
sparams.sigmay = [0,-1i;1i,0];

% % Miscellaneous simulation parameters
sparams.verbose = 1;
sparams.saveData = 1;
sparams.updateWaitbar = 3000; % How many frames to wait before updating waitbar
sparams.updateFigure = 3000; % How many frames to wait before updating figure
sparams.updateInterpPot = 40000; % Interpolate the potential at 40000 points at a time
sparams.nFigureFrames = 20; % How many frames of the simulation to save for outputted gif
sparams.nStoreDataFrames = 500; % How many times to save data

% % Stark Shift Parameters
sparams.calculateStarkShift = 0; % Set to 0 if you don't want to calcualte Stark shift, 1 if you do.
sparams.nStarkShiftFrames = 500; % How many frames of the simulation to calculate the Stark shift
sparams.v0 = 40E9; % [Hz] Reference frequency for Stark shift calculation
sparams.n2 = 2.2*1E-9*1E-9; % [m^2/V^2] Stark shift formula quadratic coefficient

% % Pulse Parameters
sparams.tcThreshold = 0.0009; % Threshold value for the tunnel coupling
sparams.voltagesToLoad{1} = 0.1;
sparams.voltagesToLoad{2} = [0.2,0.22,0.24,0.26,0.27,0.28,0.29,0.3];
sparams.voltagesToLoad{3} = [0.2,0.22,0.24,0.26,0.27,0.28,0.29,0.3];
sparams.voltagesToLoad{4} = [0.2,0.22,0.24,0.26,0.27,0.28,0.29,0.3];
sparams.voltagesToLoad{5} = 0.1;
sparams.numOfGates = length(sparams.voltagesToLoad);
% griddedInterpolant requires at least 2 grid points to interpolate, so we
% will mark which gate voltages we can actually interpolate
sparams.interpableGates = [];
for ii = 1:sparams.numOfGates
    if length(sparams.voltagesToLoad{ii}) == 1
        continue;
    end
    sparams.interpableGates = [sparams.interpableGates, ii];
end
sparams.voltagePulseBounds = [0.2,0.3;0.2,0.3;0.2,0.3];
sparams.gatesUsedInPulse = [1,2,3]; % Specificies which gates in the geometry you will actually use in the voltage pulse
sparams.nPulsePoints = length(sparams.gatesUsedInPulse)*500;
sparams.adiabaticPulseType = 'effective';

% % Adiabatic Parameters
sparams.findAdiabaticPulse = 1; % 1 if you want to use adiabatic conditions to find a pulse
% sparams.adiabaticThreshold = [[0.5,0.5]];
sparams.adiabaticThreshold = [logspace(-.25,-3,40)', logspace(-.25,-3,40)'];
sparams.hdt = 1E-15; % The h value used to find the time derivative of the wave function
sparams.timePowerBounds = [-15,-7]; % Range of 10^N time values to search in when doing optimizing the adiabatic parameter
sparams.nnIndices = 1;

% Device Geometry Parameters
sparams.dotLocs = [-140,-70,0,70,140]*1E-9;

% % Main Shuttling Simulation Parameters
sparams.dt = 5E-16; % Time between each simulation frame [sec]
sparams.sweptParameter = 'adiabicity';
% Vector of times to do the shuttling simulation over
%sparams.totalTime = 4E-9;

% % Effective Hamiltonian parameters
sparams.includeOrbital = 1;
sparams.includeExcitedOrbital = 0;
sparams.includeValley = 1;
sparams.includeSpin = 1;
sparams.includeSecondSpin = 0;
sparams.includeT2 = 0;

% % Directory Parameters
% Where to load potentials from for simulation
sparams.potDir = 'C:\Users\bbuonaco\Documents\GitHub\Simulated Potentials\FiveGate_DSize_40_GSize_30\TEMPLATE_5Gate_Dop_1.358E15_noRGrid_';
sparams.interpPotentials = 1; % Whether or not to interpolate the potentials in time domain for simulation
sparams.interpType = 'linear'; % What type of interpolation of the potentials to do in time domain
sparams.extrapType = 'linear';
for ii = 1:length(sparams.interpableGates)
    sparams.gateLabels{ii} = sprintf('V_%d',sparams.interpableGates(ii));
end
sparams.saveDir = 'C:/Users/bbuonaco/Documents/Github/Simulation Results/Electron Shuttling/';

clearvars ii1