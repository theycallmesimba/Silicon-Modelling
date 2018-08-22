%***Physical constants***%
sparams.ee = 1.602E-19;
% sparams.me = 9.11E-31*0.191;
% sparams.hbar = 6.626E-34/(2*pi);
sparams.hbar = 1;
% sparams.eps = 8.854E-12*7.8;
% sparams.me = 9.11E-31*0.067;
sparams.me = 1;
sparams.eps = 8.854E-12*12.4;

%***Input definitions***%
sparams.potFile = 'Brandon_1660mV.fig';
sparams.verbose = 1; % Will display some figures to help with debugging

%***Geometry definitions***%
% sparams.dotLocations = [-20.0,0.0; 20,0.0]*1E-9; % Each row is a dot's
% location (x,y)
a = 4;
h = sqrt(3)/2*a;
sparams.dotLocations = [0,2/3*h;-a/2,-h/3;a/2,-h/3]*aBtom;
[sparams.nDots,~] = size(sparams.dotLocations);

%***Simulation definitions***%
% How many orbitals you will consider when buiding up the localized
% harmonic orbitals.  Up to the sparams.nLocalOrbitals-1 mode will be
% considered as you always have ground state of n = 0. 
% s-orbital = 1
% p-orbital = 2
% d-orbital = 3
% etc.
sparams.nLocalOrbitals = 8;
% Now get the full number of single electron orbitals we will consider when
% we include orbitals from all dots
sparams.nSingleOrbitals = sparams.nLocalOrbitals*(sparams.nLocalOrbitals+1)/2*sparams.nDots;
% Number of grid points in x and y directions.  Initialized to 0 but filled
% in after potential is loaded
sparams.ngridx = 0;
sparams.ngridy = 0;
% Sets the maximum value for the non shfited HO we will use when we rewrite
% the basis of shifted HOs into non shifted HOs
sparams.maxNonShiftedHOsX = 25; 
sparams.maxNonShiftedHOsY = 25;
sparams.nNonShiftedHOs = sparams.maxNonShiftedHOsX*sparams.maxNonShiftedHOsY;
% Set tolerance level for norm checks
sparams.normThreshold = 1E-14;
% Set threshold for calculating the theta integrals
sparams.thetaIntThreshold = 1E-30;



