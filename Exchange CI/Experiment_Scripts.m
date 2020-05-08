
%
%*************************************************************************%
%% Convergence of Coulomb matrix elements
%% Experiment
sparams.unitsType = 'Rydberg';

maxOrgOrbs = 12;

storedHOriginB = zeros(maxOrgOrbs-1,3,3);
storedHOriginAB = zeros(maxOrgOrbs-1,3,3);
storedHLOHO = zeros(maxOrgOrbs-1,3,3);
storedHLOHOA = zeros(maxOrgOrbs-1,3,3);
storedCMEsLOHO = zeros(maxOrgOrbs-1,3^2,3^2);
storedCMEsItin = zeros(maxOrgOrbs-1,3^2,3^2);

for aa = 2:maxOrgOrbs
    fprintf(1,'*Checking convergance with orbital number %d*\n',aa);

    clearvars -except maxOrgOrbs aa stored*
    sparams.unitsType = 'Rydberg';
    
    simparams;
    sparams.maxOriginHOsX = aa; 
    sparams.maxOriginHOsY = aa;
    sparams.nOriginHOs = sparams.maxOriginHOsX*sparams.maxOriginHOsY;
    
    exCIMain2D_3bases;
    
    storedHOriginB(aa-1,:,:) = HoriginB;
    storedHOriginAB(aa-1,:,:) = HoriginAB;
    storedHLOHO(aa-1,:,:) = HLOHO;
    storedHLOHOA(aa-1,:,:) = HLOHOA;
    storedCMEsLOHO(aa-1,:,:) = CMEsHOs;
    storedCMEsItin(aa-1,:,:) = CMEsItin;
    
    fprintf(1,'*Done checking convergance with orbital number %d!*\n',aa);
end

% clearvars -except maxOrgOrbs ii stored*
% save converganceTest.mat
%% Analysis which compares calculated LOHO elements with ones from Marek's paper
load('converganceTest.mat');
startInd = 2;
endInd = maxOrgOrbs-1;

figure('Color','white','Position', [100 150 1500 800]);
%<rr|v|rr>
expectedValue = 2.939179;
subplot(2,3,1);
hold on;
set(gca,'Fontsize',14,'TicklabelInterpreter','latex');
xRange = 2:maxOrgOrbs;
inds = [1,1;5,5;9,9];
[points,~] = size(inds);
for ii = 1:points
    percentDiff = (storedCMEsLOHO(end,inds(ii,1),inds(ii,2)) - expectedValue)/expectedValue;
    plot(xRange(startInd:endInd),squeeze(storedCMEsLOHO(startInd:endInd,inds(ii,1),inds(ii,2))),...
        'Linewidth',2,'DisplayName',sprintf('%.2E',percentDiff));
end
lgd = legend('show');
lgd.FontSize = 8;
line([xRange(startInd), xRange(endInd)], [expectedValue, expectedValue],'Linestyle','--','Color','k','Linewidth',2);
title('<rr|v|rr>');
xlabel('Nx, Ny orbital states','Fontsize',18,'Interpreter','latex');
ylabel('Energy [Ry*]','Fontsize',18,'Interpreter','latex');
xlim([xRange(startInd), xRange(endInd)]);

%<rs|v|sr>
expectedValue = 0.512857;
subplot(2,3,2);
hold on;
set(gca,'Fontsize',14,'TicklabelInterpreter','latex');
xRange = 2:maxOrgOrbs;
inds = [2,2;3,3;4,4;6,6;7,7;8,8];
[points,~] = size(inds);
for ii = 1:points
    percentDiff = (storedCMEsLOHO(end,inds(ii,1),inds(ii,2)) - expectedValue)/expectedValue;
    plot(xRange(startInd:endInd),squeeze(storedCMEsLOHO(startInd:endInd,inds(ii,1),inds(ii,2))),...
        'Linewidth',2,'DisplayName',sprintf('%.2E',percentDiff));
end
lgd = legend('show');
lgd.FontSize = 8;
line([xRange(startInd), xRange(endInd)], [expectedValue, expectedValue],'Linestyle','--','Color','k','Linewidth',2);
title('<rs|v|sr>');
xlabel('Nx, Ny orbital states','Fontsize',18,'Interpreter','latex');
ylabel('Energy [Ry*]','Fontsize',18,'Interpreter','latex');
xlim([xRange(startInd), xRange(endInd)]);

%<rr|v|rs>
expectedValue = 0.004653;
subplot(2,3,3);
hold on;
set(gca,'Fontsize',14,'TicklabelInterpreter','latex');
xRange = 2:maxOrgOrbs;
inds = [1,2;1,3;1,4;1,7;2,5;3,9;4,5;6,9;7,9;8,9];
[points,~] = size(inds);
for ii = 1:points
    percentDiff = (storedCMEsLOHO(end,inds(ii,1),inds(ii,2)) - expectedValue)/expectedValue;
    plot(xRange(startInd:endInd),squeeze(storedCMEsLOHO(startInd:endInd,inds(ii,1),inds(ii,2))),...
        'Linewidth',2,'DisplayName',sprintf('%.2E',percentDiff));
end
lgd = legend('show');
lgd.FontSize = 8;
line([xRange(startInd), xRange(endInd)], [expectedValue, expectedValue],'Linestyle','--','Color','k','Linewidth',2);
title('<rr|v|rs>');
xlabel('Nx, Ny orbital states','Fontsize',18,'Interpreter','latex');
ylabel('Energy [Ry*]','Fontsize',18,'Interpreter','latex');
xlim([xRange(startInd), xRange(endInd)]);

%<rs|v|st>
expectedValue = 0.002446;
subplot(2,3,4);
hold on;
set(gca,'Fontsize',14,'TicklabelInterpreter','latex');
xRange = 2:maxOrgOrbs;
inds = [2,3;2,8;3,6;4,6;4,7;7,8];
[points,~] = size(inds);
for ii = 1:points
    percentDiff = (storedCMEsLOHO(end,inds(ii,1),inds(ii,2)) - expectedValue)/expectedValue;
    plot(xRange(startInd:endInd),squeeze(storedCMEsLOHO(startInd:endInd,inds(ii,1),inds(ii,2))),...
        'Linewidth',2,'DisplayName',sprintf('%.2E',percentDiff));
end
lgd = legend('show');
lgd.FontSize = 8;
line([xRange(startInd), xRange(endInd)], [expectedValue, expectedValue],'Linestyle','--','Color','k','Linewidth',2);
title('<rs|v|st>');
xlabel('Nx, Ny orbital states','Fontsize',18,'Interpreter','latex');
ylabel('Energy [Ry*]','Fontsize',18,'Interpreter','latex');
xlim([xRange(startInd), xRange(endInd)]);

%<rs|v|rs>
expectedValue = 0.000049;
subplot(2,3,5);
hold on;
set(gca,'Fontsize',14,'TicklabelInterpreter','latex');
xRange = 2:maxOrgOrbs;
inds = [1,5;1,9;2,4;3,7;5,9;6,8];
[points,~] = size(inds);
for ii = 1:points
    percentDiff = (storedCMEsLOHO(end,inds(ii,1),inds(ii,2)) - expectedValue)/expectedValue;
    plot(xRange(startInd:endInd),squeeze(storedCMEsLOHO(startInd:endInd,inds(ii,1),inds(ii,2))),...
        'Linewidth',2,'DisplayName',sprintf('%.2E',percentDiff));
end
lgd = legend('show');
lgd.FontSize = 8;
line([xRange(startInd), xRange(endInd)], [expectedValue, expectedValue],'Linestyle','--','Color','k','Linewidth',2);
title('<rs|v|rs>');
xlabel('Nx, Ny orbital states','Fontsize',18,'Interpreter','latex');
ylabel('Energy [Ry*]','Fontsize',18,'Interpreter','latex');
xlim([xRange(startInd), xRange(endInd)]);

%<rr|v|st>
expectedValue = 0.000019;
subplot(2,3,6);
hold on;
set(gca,'Fontsize',14,'TicklabelInterpreter','latex');
xRange = 2:maxOrgOrbs;
inds = [1,6;1,8;2,6;2,7;2,9;3,4;3,5;3,8;4,8;4,9;5,7;6,7];
[points,~] = size(inds);
for ii = 1:points
    percentDiff = (storedCMEsLOHO(end,inds(ii,1),inds(ii,2)) - expectedValue)/expectedValue;
    plot(xRange(startInd:endInd),squeeze(storedCMEsLOHO(startInd:endInd,inds(ii,1),inds(ii,2))),...
        'Linewidth',2,'DisplayName',sprintf('%.2E',percentDiff));
end
lgd = legend('show');
lgd.FontSize = 8;
line([xRange(startInd), xRange(endInd)], [expectedValue, expectedValue],'Linestyle','--','Color','k','Linewidth',2);
title('<rr|v|st>');
xlabel('Nx, Ny orbital states','Fontsize',18,'Interpreter','latex');
ylabel('Energy [Ry*]','Fontsize',18,'Interpreter','latex');
xlim([xRange(startInd), xRange(endInd)]);
%%
%
%*************************************************************************%
%% Time comparison of optimized CME calculator and non-optimized CME calculator
%% Experiment
startOrbs = 6;
endOrbs = 6;
timings = zeros(2,length(startOrbs:endOrbs)); %ind 1 = non-optimized, ind 2 = optimized
% timings = [];
kk = 0;
for ii = startOrbs:endOrbs
    fprintf('$**********************$\n');
    fprintf('Number of orbitals = %d\n', ii);
    fprintf('$**********************$\n');
    
    [~, kk] = size(timings);
    kk = kk + 1;
    
    sparams.maxOriginHOsX = ii;
    sparams.maxOriginHOsY = ii;
    sparams.nOriginHOs = sparams.maxOriginHOsX*sparams.maxOriginHOsY;
    sparams.originHOs = createOriginHOs(sparams, gparams, 1.3749, 0);
    
%     fprintf(1,'Running non-optimized code...\n');
%     tic;
%     temp0 = constructOriginCMEsLib(sparams, 1);
%     timings(1,kk) = toc;
%     fprintf(1,'Done! Took %d seconds\n',round(timings(1,kk)));
%     
    fprintf(1,'Running lookup code...\n');
    tic;
    temp1 = constructOriginCMEsLibOPTIMIZED(sparams, 1);
    timings(2,kk) = toc;
    fprintf(1,'Done! Took %d seconds\n',round(timings(2,kk)));
    
%     fprintf(1,'Running GPU code...\n');
%     tic;
%     temp2 = constructOriginCMEsLibGPU(sparams, 1);
%     timings(3,kk) = toc;
%     fprintf(1,'Done! Took %d seconds\n',round(timings(3,kk)));
    
    fprintf('Speedup (GPU/REG) = %.2E\n', timings(3,kk)/timings(1,kk));
    fprintf('Speedup (LKUP/REG) = %.2E\n', timings(2,kk)/timings(1,kk));
    fprintf('Speedup (GPU/LKUP) = %.2E\n', timings(3,kk)/timings(2,kk));
%     temp0 = full(temp0);
end
%% Analysis of experiment
figure;
yyaxis left
semilogy(1:11, timings);
ylabel('Time [sec]','Interpreter','latex','Fontsize',20);
yyaxis right
plot(1:11, timings);
ylabel('Time [sec]','Interpreter','latex','Fontsize',20);
set(gca,'TickLabelInterpreter','latex','Fontsize',14);
xlabel('$N_x$, $N_y$ harmonic orbitals','Fontsize',20,'Interpreter','latex');
[~,maxOrb] = size(timings);
xlim([1,maxOrb]);
%%
%
%*************************************************************************%
%% Energy convergance with omega optimization and no omega optimization
%% Experiment using intermediate LOHO basis
nOrbsStart = 3;
nOrbsEnd = 14;
manyBodyEns = zeros(nOrbsEnd-nOrbsStart+1,3);
manyBodyEnsOptOmega = zeros(nOrbsEnd-nOrbsStart+1,3);
for jj = [0,1]
    for ii = nOrbsStart:nOrbsEnd
        Initialize_Potential;
%         sparams.maxOriginHOsX = ii;
%         sparams.maxOriginHOsY = ii;
%         sparams.nItinerantOrbitals = 3;
        sparams.maxOriginHOsX = 11;
        sparams.maxOriginHOsY = 11;
        sparams.nItinerantOrbitals = ii;
        sparams.nOriginHOs = sparams.maxOriginHOsX*sparams.maxOriginHOsY;
        [eVectors, ens] = calculateManyBodySpectra_3Bases(sparams, gparams, jj, 0);
%         [eVectors, ens] = calculateManyBodySpectra_2Bases(sparams, gparams, jj, 0);
        diag(ens)
        if jj == 0
            manyBodyEns(ii,:) = diag(ens(1:3,1:3));
        elseif jj == 1
            manyBodyEnsOptOmega(ii,:) = diag(ens(1:3,1:3));
        end
        clearvars sparams;
    end
end
%% Experiment analysis
figure;
hold on;
hG1Lines = plot(nOrbsStart:nOrbsEnd,manyBodyEns(nOrbsStart:nOrbsEnd,:),'k','linewidth',1.5);
hG2Lines = plot(nOrbsStart:nOrbsEnd,manyBodyEnsOptOmega(nOrbsStart:nOrbsEnd,:),'r','linewidth',1.5);
set(gca,'TickLabelInterpreter','latex','Fontsize',14);
xlabel('Nx, Ny origin HOs','Interpreter','latex','Fontsize',20);
xlabel('N canonical orbitals','Interpreter','latex','Fontsize',20);
ylabel('Energy [Ry*]','Interpreter','latex','Fontsize',20);
xlim([nOrbsStart,nOrbsEnd]);
% xlim([9,nOrbsEnd]);
% ylim([-24.5,-23.9]);
% title('Energy convergence: 10 Canons','Interpreter','latex','Fontsize',20);
title('Energy convergence','Interpreter','latex','Fontsize',20);
hG1Group = hggroup;
hG2Group = hggroup;
set(hG1Lines,'Parent',hG1Group)
set(hG2Lines,'Parent',hG2Group)
% Include these hggroups in the legend:
set(get(get(hG1Group,'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','on'); 
set(get(get(hG2Group,'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','on'); 
legend({'without $\omega$ optimization','with $\omega$ optimization'},'Location','northeast','Interpreter','latex')
legend('boxoff')
% export_fig 'EConvWRTw_Canon-10' -m5
%%
%
%*************************************************************************%


sparams.maxOriginHOsX = 7;
sparams.maxOriginHOsY = 7;
sparams.nOriginHOs = sparams.maxOriginHOsX*sparams.maxOriginHOsY;
sparams.originHOs = createOriginHOs(sparams, gparams, 1.3749, 0);

fprintf(1,'Running GPU code...\n');
tic;
temp0 = full(constructOriginCMEsLibGPU(sparams, 1));
timings = toc;
fprintf(1,'Done! Took %d seconds\n',round(timings));






