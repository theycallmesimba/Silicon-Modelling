vArray = [0.8,0.798416157705339,0.6]';
currPot = sparams.P2DEGInterpolant(getInterpolantArgument(vArray,xx));
currPot = squeezeFast(sparams.numOfGates,currPot)';
tcs = calculateTunnelCoupling(sparams, xx, currPot)/sparams.ee

[wfs,ens] = solve1DSingleElectronSE(sparams, 2, xx, currPot);
tcsOrig = (ens(2,2) - ens(1,1))/2/sparams.ee

figure;
hold on;
yyaxis left
plot(xx,currPot);
yyaxis right
plot(xx,wfs.^2);

%%

% T2Sweep = logspace(-10,-8,20)/1E-9;
tcSweep = logspace(-7,-4,10)/1E-6;
detuningSweep = logspace(0,2,8);

figure;
[T2XX,tcYY] = meshgrid(detuningSweep,tcSweep);
s = surf(T2XX,tcYY,real(purityZZ));
set(gca,'Fontsize',16);
set(gca,'YScale','log','XScale','log');
% shading interp
cbar = colorbar;
cbar.Label.String = 'Purity';
cbar.Label.Interpreter = 'latex';
cbar.Label.FontSize = 20;
view(2);
xlabel('$\max\epsilon$ [ns]','Interpreter','Latex','Fontsize',20);
ylabel('$t_c \,\, {\rm [\mu eV]}$','Interpreter','Latex','Fontsize',20);

% figure;
% s = surf(T2XX,tcYY,T2Overtc);
% set(gca,'Fontsize',16);
% set(gca,'XScale','log','YScale','log');
% shading interp
% cbar = colorbar;
% cbar.Label.String = '$T_2/t_c$';
% cbar.Label.Interpreter = 'latex';
% cbar.Label.FontSize = 20;
% view(2);
% xlabel('T2 [ns]','Interpreter','Latex','Fontsize',20);
% ylabel('$t_c \,\, {\rm [\mu eV]}$','Interpreter','Latex','Fontsize',20);