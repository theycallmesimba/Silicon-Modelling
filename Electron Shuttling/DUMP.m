shuttleParameterFile;
currFName = 'V1_0.100_V2_0.200_V3_0.300_V4_0.300_V5_0.100/output/potential';
[xx, zz, pot2D_XZ] = loadPotentialFile(sparams,[sparams.potDir currFName],...
    1,'nextnano','x',0);

figure;
hold on;
[~,twoDEGindZ] = min(abs(zz - (-1*1E-9)));
plot(xx,pot2D_XZ(twoDEGindZ,:));
[~,twoDEGindZ] = min(abs(zz - (-3*1E-9)));
plot(xx,pot2D_XZ(twoDEGindZ,:));