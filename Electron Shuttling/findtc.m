shuttleParameterFile;

% gapSizes = [7,9,11,13,15]; % [nm]
% gapSizes = [10,15,20];
gapSizes = 7;
dotSizes = 150; %[nm]

potDirectoryLocation = sparams.potDir;

nPts = 2;
tcs = zeros(length(gapSizes),nPts);
pTime = zeros(length(gapSizes),nPts);
% Record what the actual bounds we found are and the actual tunnel
% couplings we can find in the range.
vBoundsActual = zeros(length(gapSizes),2);
tcBoundsActual = zeros(length(gapSizes),2);
    
for jj = 1:length(gapSizes)
    fprintf(1,'Analyzing gapsize %d potentials...\n',gapSizes(jj));
    sparams.potDir = [potDirectoryLocation, '150nm ' num2str(gapSizes(jj)) 'nm-gap\'];

%     sparams.voltagesToLoad{1} = 0.18:0.02:0.22;
%     sparams.voltagesToLoad{2} = 0.18:0.02:0.22;
%     sparams.voltagesToLoad{3} = sparams.voltagesToLoad{1} - 0.2;
%     sparams.voltagesToLoad{1} = 0.41:0.05:0.61;
%     sparams.voltagesToLoad{2} = 0.41:0.05:0.61;
%     sparams.voltagesToLoad{3} = 0.41:0.05:0.61;
    sparams.voltagesToLoad{1} = 0.15:0.05:0.55;
    sparams.voltagesToLoad{2} = 0.15:0.05:0.55;
    sparams.voltagesToLoad{3} = 0.15:0.05:0.55;

    fprintf(1,'Loading potentials...\n');
    [sparams,xx,zz] = loadPotentials(sparams);
    sparams.nxGrid = length(xx);
    sparams.nzGrid = length(zz);
    sparams.dx = xx(2) - xx(1);
    sparams.dz = zz(2) - zz(1);

    % Find which index corresponds to where the 2DEG should be
    [~,sparams.twoDEGindZ] = min(abs(zz - (-0.5*1E-9)));
    for ii = 1:length(sparams.potentials)
        sparams.potentials(ii).pot2DEG = sparams.potentials(ii).pot2D(sparams.twoDEGindZ,:);
    end

    % Now we want to make the potential interpolant object (both 2D and 2DEG)
    sparams.interpType = 'spline';
    sparams.extrapType = 'spline';
    sparams = makePotentialsInterpolants(sparams,xx,zz);

    dotLocs = [-dotSizes - gapSizes(jj),0,dotSizes + gapSizes(jj)]*1E-9;
    
    vProbe = 0.55;
    vVec = [vProbe,vProbe,vProbe-0.2];
    [vgSwept, tcProbe] = findResonantTunnelCoupling(sparams, xx, vVec, 2, dotLocs, 1);
    
    fprintf(1,'Found tc bounds %0.4E [eV] with V %0.6E [V]\n',tcProbe/sparams.ee,vgSwept);
end

%%
vProbe = 0.55;
vSweep = 0.5100204;
vVec = [vProbe,vSweep,vProbe-0.2];
plotPotentialAndGroundWF( sparams, vVec, xx );

%%
% v1 = 0.5100204;
% v2 = 0.5100204;
v1 = 1.5;
v2 = 1.7;
v3 = 1.5;
gVolts = [v1, v2, v3];
plotPotentialAndGroundWF( sparams, gVolts, xx );
tempPot = sparams.P2DEGInterpolant([num2cell(gVolts),mat2cell(xx,1,length(xx))]);
tempPot = squeezeFast(sparams.numOfGates,tempPot);
[~, ens] = solve1DSingleElectronSE(sparams,2,xx,tempPot); 

fprintf(1,'Orbital spacing DeltaE = %.8E [eV]\n',abs(ens(2,2) - ens(1,1))/sparams.ee);
%%
xx = [40, 80, 150];
yy = [3.3E-3, 9.3E-4, 3.7E-4];
plot(xx,yy);
