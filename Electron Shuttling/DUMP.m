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

%     sparams.voltagesToLoad{1} = 0.2:0.25:2.0;
%     sparams.voltagesToLoad{2} = 0.2:0.25:2.0;
%     sparams.voltagesToLoad{3} = sparams.voltagesToLoad{1} - 0.2;
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
    
    tcBounds = [1E-5,1E-4]*sparams.ee;
    vSearchBnds = [min(sparams.voltagesToLoad), max(sparams.voltagesToLoad)];

    vBounds = findVoltageBoundsForTunnelCouplingBounds( sparams, xx, tcBounds, vSearchBnds, dotLocs);  
    
    vBoundsActual(jj,:) = vBounds;
    vVec = [vBounds(1),vBounds(1),vBounds(1)-0.2];
    [~, tcBoundsActual(jj,1)] = findResonantTunnelCoupling(sparams, xx, vVec, 2, dotLocs, 1);
    vVec = [vBounds(2),vBounds(2),vBounds(2)-0.2];
    [~, tcBoundsActual(jj,2)] = findResonantTunnelCoupling(sparams, xx, vVec, 2, dotLocs, 1);
    
    fprintf(1,'Found tc bounds [%.4E,%.4E] [eV] with corresponding voltages [%.4E,%.4E] [V]\n',...
        tcBoundsActual(jj,1)/sparams.ee,tcBoundsActual(jj,2)/sparams.ee,...
        vBoundsActual(jj,1),vBoundsActual(jj,2));
    
    V = linspace(vBounds(1),vBounds(2),nPts);
    V = [];

%     figure;
%     title(sprintf('Gap size = %d',gapSizes(jj)));  
%     hold on;

    for ii = 1:length(V)
        vVec = [V(ii),V(ii),V(ii)-0.2];
        
        [vGTargMax, tcs(jj,ii)] = findResonantTunnelCoupling(sparams, xx, vVec, 2, dotLocs, 1);
%         fprintf(1,'V = %0.4E [V], tc = %0.6E [eV]\n',V(ii),tcs(jj,ii)/sparams.ee);
     
%         vVec(2) = vGTargMax;
%         vVecCell = [num2cell(vVec),xx];
%         currPot = squeezeFast(sparams.numOfGates,sparams.P2DEGInterpolant(vVecCell));
%         plot(xx,currPot/sparams.ee);
        
        
%       [~, ~, pTime(jj,ii)] = getVoltagePulseAdiabatic(sparams,xx,[0.001,0.001],vGTargMax,dotLocs);
%       fprintf(1,'tc = %0.6E [eV], pulse time = %0.6E [s]\n',tcs(jj,ii)/sparams.ee,pTime(jj,ii));
    end
    
%     figure;
%     plot(V,tcs);
end
%%
figure;
hold on;
for ii = 1:length(gapSizes)
    plot(tcs(ii,:)/sparams.ee,(dotSizes + gapSizes(ii))*1E-9./pTime(ii,:),'Linewidth',2);
    legVec{ii} = ['DotSize-' num2str(dotSizes) 'nm GapSize-' num2str(gapSizes(ii)) 'nm'];
end
set(gca,'Fontsize',15);
xlabel('Tunnel Coupling [eV]','Interpreter','Latex','Fontsize',25);
ylabel('Electron Velocity [nm/ns]','Interpreter','Latex','Fontsize',25);
legend(legVec)
 figure;;
 yyaxis left
 plot(V,tcs/sparams.ee); 
 yyaxis right
 plot(V,diff/sparams.ee);
 
 figure;
 plot(tcs/sparams.ee, pTime);

% plotPotentialAndGroundWF( sparams, [2,2.09,1.8], xx );
%  for ii = 1: length(V);
%     diff(ii) = findPeakDiff(V(ii),sparams, xx, ...
%       2,[V(ii),V(ii)+0.2,V(ii)-0.2]);
%  end
%  figure;
%  plot(V,diff);
% for ii =1:length(V);
%     
% [wfs, ens ] = solve1DSingleElectronSE ( sparams, 2, X, V(ii) );
% end
% x1 = [1.57E-23	1.48E-23	1.39E-23	1.31E-23	1.22E-23	1.14E-23	1.07E-23	9.94E-24	9.25E-24	8.60E-24	7.98E-24  7.39E-24	6.84E-24	6.33E-24	5.85E-24	5.40E-24	4.98E-24	4.59E-24	4.23E-24	3.89E-24	3.58E-24	3.29E-24	3.02E-24	2.77E-24	2.55E-24	2.34E-24	2.14E-24	1.96E-24	1.80E-24	1.65E-24];
% y1 = [3.88E-09	4.17E-09	4.50E-09	4.85E-09	5.24E-09	5.66E-09	6.15E-09	6.64E-09	7.20E-09	7.85E-09	8.59E-09  9.24E-09	1.01E-08	1.09E-08	1.21E-08	1.30E-08	1.41E-08	1.58E-08	1.72E-08	1.83E-08	2.05E-08	2.16E-08	2.24E-08	2.03E-08	2.10E-08	2.03E-08	2.06E-08	2.02E-08	2.02E-08	1.96E-08];
% x3 = [5.81E-24	5.35E-24	4.92E-24	4.52E-24	4.14E-24	3.79E-24	3.47E-24	3.17E-24	2.89E-24	2.64E-24	2.40E-24  2.19E-24	1.99E-24];
% y3 = [1.29E-08	1.40E-08	1.54E-08	1.69E-08	1.87E-08	2.06E-08	2.26E-08	2.48E-08	2.54E-08	2.64E-08	2.65E-08  2.60E-08	2.64E-08];
% figure;
% hold on;
% set(gca,'Fontsize',16);
% plot(x1/sparams.ee,(7 + 60)./y1*1E-9,'Linewidth',1.5);
% plot (tcs40nmDot10nmGap/sparams.ee,pTime40nmDot10nmGap);
% plot(x3/sparams.ee,(7 + 80)./y3*1E-9,'Linewidth',1.5);
% plot(tt/sparams.ee,(10 + 80)./pp*1E-9,'Linewidth',1.5);
% plot(D60nmG10nm/sparams.ee,(10 + 60)./D60nmG10nmpTime*1E-9,'Linewidth',1.5);
% ylabel('Voltage Pulse Time [s]','Interpreter','Latex','Fontsize',25);
% ylabel('Electron Velocity [nm/ns]','Interpreter','Latex','Fontsize',25);
% xlabel('Tunnel Coupling  [eV]','Interpreter','Latex','Fontsize',25);
% legend('60nm DOT 7nm GAP','80nm DOT 7nm GAP','80nm DOT 10nm GAP','60nm DOT 10nm GAP','Location','northwest');
% 
% save('TUNNELDATA','x1','y1','x3','y3','tt','pp','D60nmG10nm','D60nmG10nmpTime');