%***Physical constants***%
sparams.ee = 1.602E-19;
sparams.me = 9.11E-31*0.191;
sparams.hbar = 6.626E-34/(2*pi);
sparams.eps = 8.854E-12*7.8;

%***Input definitions***%
sparams.potFile = 'Brandon_1660mV.fig';
sparams.verbose = 1; % Will display some figures to help with debugging

%***Geometry definitions***%
sparams.dotLocations = [-20.0,0.0; 20,0.0]*1E-9; % Each row is a dot's location x,y
[sparams.nDots,~] = size(sparams.dotLocations);

%***Simulation definitions***%
% How many orbitals you will consider when buiding up the localized
% harmonic orbitals.  Up to the sparams.nLocalOrbitals-1 mode will be
% considered as you always have ground state of n = 0
sparams.nLocalOrbitals = 12;
% This parameter dictates how many combinations of the x and y localized
% HOs you use as your basis set.  maxLocal2DOrbital = 2 would imply the
% following (n,m) pairs are taken as your basis set {(0,0), (0,1), (0,2),
% (1,0), (1,1), (2,0)}.  Note that (1,2), (2,1), and (2,2) are omitted.
sparams.maxLocal2DOrbital = 5;
% Total number of single particle orbitals we will have found after
% solving for the 1 electron wave function and after accounting for
% maxLocal2DOrbital
an = 1;
for ii = 1:sparams.maxLocal2DOrbital
    an = an + ii + 1;
end
sparams.nSingleOrbitals = an*sparams.nDots;
clearvars an ii;
% Number of grid points in x and y directions.  Initialized to 0 but filled
% in after potential is loaded
sparams.ngridx = 0;
sparams.ngridy = 0;
% Sets the maximum value for the non shfited HO we will use when we rewrite
% the basis of shifted HOs into non shifted HOs
sparams.maxNonShiftedHOsX = 3; 
sparams.maxNonShiftedHOsY = 3;
sparams.nNonShiftedHOs = ((sparams.maxNonShiftedHOsX+1)*(sparams.maxNonShiftedHOsY+1));
% Set tolerance level for norm checks
sparams.normThreshold = 1E-14;
% Set threshold for calculating the theta integrals
sparams.thetaIntThreshold = 1E-30;



