%***Physical constants***%
sparams.unitsType = 'Rydberg';
sparams.materialSystem = 'GaAs';
if strcmp(sparams.unitsType,'SI')
    sparams.hbar = 6.626E-34/(2*pi);
    if strcmp(sparams.materialSystem,'Silicon')
        sparams.eps = 8.854E-12*7.8;
        sparams.me = 9.11E-31*0.191;
    elseif strcmp(sparams.materialSystem,'GaAs')
        sparams.eps = 8.854E-12*12.4;
        sparams.me = 9.11E-31*0.067;      
    end
elseif strcmp(sparams.unitsType,'Rydberg')
    sparams.hbar = 1;
    sparams.me = 1;
    if strcmp(sparams.materialSystem,'Silicon')
        sparams.eps = 8.854E-12*7.8;
    elseif strcmp(sparams.materialSystem,'GaAs')
        sparams.eps = 8.854E-12*12.4;
    end
end
sparams.ee = 1.602E-19;

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
sparams.nLocalOrbitals = 1;
% Now get the full number of single electron orbitals we will consider when
% we include orbitals from all dots
sparams.nSingleOrbitals = sparams.nLocalOrbitals*(sparams.nLocalOrbitals+1)/2*sparams.nDots;
% Number of grid points in x and y directions.  Initialized to 0 but filled
% in after potential is loaded
sparams.ngridx = 0;
sparams.ngridy = 0;
% Sets the maximum value for the non shfited HO we will use when we rewrite
% the basis of shifted HOs into non shifted HOs
sparams.maxOriginHOsX = 5; %10; 
sparams.maxOriginHOsY = 5; %10;
sparams.nOriginHOs = sparams.maxOriginHOsX*sparams.maxOriginHOsY;
% Set tolerance level for norm checks
sparams.normThreshold = 1E-14;
% Set threshold for calculating the theta integrals
sparams.thetaIntThreshold = 1E-30;



