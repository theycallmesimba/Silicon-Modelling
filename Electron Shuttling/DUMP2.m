% Load all the parameters for the simulation
clear sparams xx vv pp;
shuttleParameterFile;
%    '50nm dot 0.500-2.00V','60nm dot 0.500-2.00V','70nm dot 0.500-2.00V',... 
%     '80nm dot 0.500-2.00V',
potDirs = {'90nm dot 0.500-2.00V','100nm dot 0.500-2.00V',...
    '110nm dot 0.500-2.00V','120nm dot 0.500-2.00V'};
dSize = 90:10:120; 
range = [0.005,0.005,0.004,0.004];

% potDirs = {'40nm dot potentials','60nm dot potentials','80nm dot potentials'};
% dSize = 40:20:80;

V = 0.5:0.15:2.0;
tcs = zeros(length(V),length(dSize));
diff = zeros(length(V),length(dSize));
for bb = 1:length(dSize)
    sparams.potDir = ['C:\Users\bbuonaco\Documents\MATLAB\Simulated Potentials\changing dot potentials\'...
        potDirs{bb} '/'];
    
    fprintf(1,'(%d/%d) Loading potentials for dotsize %d...\n',bb,length(dSize),dSize(bb));
    [sparams,xx,zz] = loadPotentials(sparams);
    
    sparams.nxGrid = length(xx);
    sparams.nzGrid = length(zz);
    sparams.dx = xx(2) - xx(1);
    sparams.dz = zz(2) - zz(1);
    sparams.dp = 2*pi*sparams.hbar/(sparams.dx*sparams.nxGrid);
    pp = ((-sparams.nxGrid/2):1:(sparams.nxGrid/2 - 1))*sparams.dp;

    % Find which index corresponds to where the 2DEG should be
    [~,sparams.twoDEGindZ] = min(abs(zz - (-0.5*1E-9)));
    for aa = 1:length(sparams.potentials)
        sparams.potentials(aa).pot2DEG = sparams.potentials(aa).pot2D(sparams.twoDEGindZ,:);
    end

    % Now we want to make the potential interpolant object (both 2D and 2DEG)
    % sparams = makePotentialsInterpolants(sparams,xx,zz);
    sparams = makePotentialsInterpolants(sparams,xx,zz);

    for aa = 1:length(V)
        [vGTargMax, tcs(aa,bb)] = findTunnelCouplingMax(sparams, xx, [V(aa),V(aa),V(aa)-0.2],V(aa) - range(bb),...
            V(aa) + range(bb), 2);
        diff(aa,bb) = findPeakDiff(V(aa),sparams, xx, 2,[V(aa),vGTargMax,V(aa)-0.2]); 
    end
end

[DSIZE,VV] = meshgrid(dSize,V);
figure;
title('Tc')
s = surf(DSIZE,VV,tcs/sparams.ee);
set(s,'Edgecolor','none');
colorbar;

figure;
title('DIFF')
s = surf(DSIZE,VV,diff/sparams.ee);
set(s,'Edgecolor','none');
colorbar

