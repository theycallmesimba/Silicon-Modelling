shuttleParameterFile;
% folder = 'C:\Users\bbuonaco\Documents\GitHub\Simulated Potentials\Proper 2 gate simulation';
folder = 'C:\Users\bbuonaco\Documents\GitHub\Simulated Potentials\Proper 4 gate simulation';

% Vg = linspace(0,1.5,31);
gap = [30,35,40,45,50];
gap = 40;
gwidth = 40;
% grid = [1,2,3,4,5,6,7,8,9];
doping = '1.358E15';
figure;
hold on;
% sweepVar = Vg;
sweepVar = gap;
% sweepVar = grid;
for ii = 1:length(sweepVar)
%     if Vg(ii) == 0
%         exponent = 0;
%         base = 0;
%     else
%         exponent = floor(log10(Vg(ii)));
%         base = Vg(ii)*10^(-exponent);
%     end
%     currFoldName = [folder '\' sprintf('VL1_VR1_Sweep_V_L1_%.1fE%d',base,exponent)];
%     tempfold = ['\' sprintf('TEMPLATE_2Gate_V_L1_%.1fE%d',base,exponent)];
%     fileToLoad = [currFoldName tempfold '\output\potential'];
%     currFoldName = [folder '\' sprintf('TEMPLATE_2Gate_G_GAP_%d_',sweepVar(ii)) doping];
%     currFoldName = [folder '\TEMPLATE_4Gate_G_GAP_35_' doping sprintf('_coarserGrid%d',ii)];
%     currFoldName = [folder '\TEMPLATE_4Gate_G_GAP_35_' doping sprintf('_coarserRepeatGrid%d',sweepVar(ii))];
    currFoldName = [folder '\TEMPLATE_4Gate_Dop_' doping '_RGrid_G_WIDTHX_' num2str(gwidth) '_G_GAP_' num2str(sweepVar(ii))];
    fileToLoad = [currFoldName '\output\potential'];

    [xx, zz, pot2D] = loadPotentialFile(fileToLoad, 1, 'nextnano','x');
    pot2D = -sparams.ee*pot2D;
    
    % Find which index corresponds to where the 2DEG should be
    twoDEGindZ = getClosestArrayIndex(-0.6*1E-9, zz);
    pot2DEG = pot2D(twoDEGindZ, :);

    plot(xx/1E-9, pot2DEG/sparams.ee, 'DisplayName',sprintf('%d',ii));
    xlim([min(xx),max(xx)]/1E-9);
    
    [~, ens] = solve1DSingleElectronSE(sparams,2,xx,pot2DEG);
    tc = (ens(2,2) - ens(1,1))/sparams.ee/2;
    fprintf(1,'Sweep variable %d gives tc = %E\n',sweepVar(ii),tc)
end
legend;

[psi0, ens] = solve1DSingleElectronSE(sparams,3,xx,pot2DEG);
line([min(xx),max(xx)]/1E-9,[ens(1,1),ens(1,1)]/sparams.ee);
line([min(xx),max(xx)]/1E-9,[ens(2,2),ens(2,2)]/sparams.ee);
% line([min(xx),max(xx)]/1E-9,[ens(3,3),ens(3,3)]/sparams.ee);

% tc = (ens(2,2) - ens(1,1))/sparams.ee;
figure;
plot(xx/1E-9,psi0.^2);






