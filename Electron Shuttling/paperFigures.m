% Five gate simulation pictures
% Figure 1
directory = 'C:\Users\bbuonaco\Desktop\';
currFig = 'V1_0.8_V2_0.6_V3_0.6_V4_0.6_V5_0.6.fig';
fig = open([directory currFig]);
axes = get(fig,'Children');
data = get(axes,'Children');
xdata = get(data,'XData');
ydata = get(data,'YData');
% zdata = get(data,'ZData');
delete(fig);

buf = 150;
ee = 1.602E-19;
fontsize = 30;

desiredGridX = 2^(nextpow2(length(xdata)));
xxq = linspace(min(xdata),max(xdata),desiredGridX);
vvq = interp1(xdata,ydata,xxq);

tempPot = vvq(buf:(end-buf));
tempXX = xxq(buf:(end-buf));

[tempWF, ~] = solve1DSingleElectronSE(sparams,1,tempXX,tempPot); 
tempXX = tempXX/1E-9;
figure('pos',[0 0 1200 500])
set(gca,'Fontsize',18);
xlabel('Position [nm]','Interpreter','Latex','Fontsize',fontsize,'Fontweight','bold');
xlim([min(tempXX),max(tempXX)]);

yyaxis left
plot(tempXX,tempPot/ee,'Linewidth',2.5);
ylabel('Potential [V]','Interpreter','Latex','Fontsize',fontsize,'Fontweight','bold');

yyaxis right
plot(tempXX,abs(tempWF).^2/norm(abs(tempWF).^2),'Linewidth',2.5);
ylabel('Probability','Interpreter','Latex','Fontsize',fontsize,'Fontweight','bold');
%%
% Figure 2
currFig = 'V1_0.6_V2_0.6_V3_0.6_V4_0.6_V5_0.8.fig';
fig = open([directory currFig]);
axes = get(fig,'Children');
data = get(axes,'Children');
xdata = get(data,'XData');
ydata = get(data,'YData');
% zdata = get(data,'ZData');
delete(fig);

tempPot = ydata(buf:(end-buf));
tempXX = xdata(buf:(end-buf));
[tempWF, ~] = solve1DSingleElectronSE(sparams,1,tempXX,tempPot); 
tempXX = tempXX/1E-9;
figure('pos',[0 0 1200 500])
set(gca,'Fontsize',18);
xlabel('Position [nm]','Interpreter','Latex','Fontsize',fontsize,'Fontweight','bold');
xlim([min(tempXX),max(tempXX)]);

yyaxis left
plot(tempXX,tempPot/sparams.ee,'Linewidth',2.5);
ylabel('Potential [V]','Interpreter','Latex','Fontsize',fontsize,'Fontweight','bold');

yyaxis right
plot(tempXX,abs(tempWF).^2/norm(abs(tempWF).^2),'Linewidth',2.5);
ylabel('Probability','Interpreter','Latex','Fontsize',fontsize,'Fontweight','bold');
%%
% Figure 3
directory = 'C:\Users\bbuonaco\Desktop\5 Gate 40nm CSV Files\';

g2 = [7.5, 7.6, 7.7, 7.9, 8.0, 8.1];
g3 = g2;
nn = 1;
for ii = 1:length(g2)
    for jj = 1:length(g3)
        currString = sprintf('V1_0.6_V2_0.6_V3_0.6_V4_0.6_V5_0.6_V_L1_%0.1f_V_C_%0.1f.csv',g2(ii),g3(jj));
        
        data = dlmread([directory currString]);
        [rows,cols] = size(data);
        zdata = data(2:rows,2:cols);
        
        % Then we are on the first potential
        if nn == 1
            xdata = data(1,2:cols);
            ydata = data(2:rows,1);

            % The data may not be uniform in sampling, so we need to fix that for
            % the fourier transforms in the main code for speed up.
            % Find next highest power of 2 to the length of xx
            desiredGridX = 2^(nextpow2(length(xdata)));
            desiredGridZ = round(2.5*length(ydata));

            % Make linearly spaced grid of points to interpolate the
            % potential at
            xxq = linspace(min(xdata),max(xdata),desiredGridX);
            zzq = linspace(min(ydata),max(ydata),desiredGridZ);

            [XX,ZZ] = meshgrid(xdata,ydata);
            [XXq,ZZq] = meshgrid(xxq,zzq);
%                 pots = zeros(1,desiredGridZ,desiredGridX);
        end
        potentials(nn).pot2D = interp2(XX,ZZ,zdata,XXq,ZZq);
        potentials(nn).pot2D = -potentials(nn).pot2D*ee; % Convert to J
        potentials(nn).gateValues = [g2(ii) g3(jj)];
        nn = nn + 1;
    end
end
%%
% Find which index corresponds to where the 2DEG should be
[~,twoDEGindZ] = min(abs(zzq - (-0.5*1E-9)));
for ii = 1:length(potentials)
    potentials(ii).pot2DEG = potentials(ii).pot2D(twoDEGindZ,:);
end

numOfGates = 2;
    
% Get the min and max value for each gate
allGateValues = reshape([potentials.gateValues],numOfGates,[])';

g2GVec = unique(allGateValues(:,1));
g3GVec = unique(allGateValues(:,2));

nzGrid = length(zzq);
nxGrid = length(xxq);

% Create the potentials interpolant
P2DInterpolant = griddedInterpolant({g2GVec,g3GVec,zzq,xxq},...
    permute(reshape([potentials.pot2D],...
    [nzGrid,nxGrid,length(g3GVec),length(g2GVec)]),[4,3,1,2]));

P2DEGInterpolant = griddedInterpolant({g2GVec,g3GVec,xxq},...
    permute(reshape([potentials.pot2DEG],...
    [nxGrid,length(g3GVec),length(g2GVec)]),[3,2,1]));

testPot = squeeze(P2DEGInterpolant({8.0,7.9978,xxq}))';

tempPot = testPot(buf:(end-buf));
tempXX = xxq(buf:(end-buf));
tempXX = tempXX*1E-9;

[tempWF, ~] = solve1DSingleElectronSE(sparams,1,tempXX,tempPot);
tempXX = tempXX/1E-9;
figure('pos',[0 0 1200 500])
set(gca,'Fontsize',18);
xlabel('Position [nm]','Interpreter','Latex','Fontsize',fontsize,'Fontweight','bold');
xlim([min(tempXX),max(tempXX)]);

yyaxis left
plot(tempXX,tempPot/sparams.ee,'Linewidth',2.5);
ylabel('Potential [V]','Interpreter','Latex','Fontsize',fontsize,'Fontweight','bold');

yyaxis right
plot(tempXX,abs(tempWF).^2/norm(abs(tempWF).^2),'Linewidth',2.5);
ylabel('Probability','Interpreter','Latex','Fontsize',fontsize,'Fontweight','bold');

%%
% Now for the big 3D model picture
file = 'C:\Users\bbuonaco\Desktop\potential.dat';
data = dlmread(file);
data = reshape(data,[377,69,42]);

file = 'C:\Users\bbuonaco\Desktop\potential.coord';
coords = dlmread(file);
xx = coords(1:377);
yy = coords(378:446);
zz = coords(447:end);
bufX = 80;
bufY = 1;
[XX,YY] = meshgrid(xx(bufX:(end-bufX)),yy(bufY:(end-bufY)));

% Find which index corresponds to where the 2DEG should be
[~,twoDEGindZ] = min(abs(zz - (-0.5)));
twoDPot2DEG = -squeeze(data(bufX:(end-bufX),bufY:(end-bufY),twoDEGindZ))';
twoDPot2DEG = twoDPot2DEG - max(max(twoDPot2DEG));
twoDPot2DEG = twoDPot2DEG/abs(min(min(twoDPot2DEG)));
twoDPot2DEG = twoDPot2DEG*120;

s = surf(XX,YY,twoDPot2DEG);
% axis vis3d
% h1 = light;
% h2 = light;
% h3 = light;
% s.FaceLighting = 'gouraud';
% s.AmbientStrength = 0.4;
% s.DiffuseStrength = 0.8;
% s.SpecularStrength = 0.9;
% s.SpecularExponent = 25;
% s.BackFaceLighting = 'unlit';
% camlight('left')

hold on;
set(s,'EdgeColor','none');
xlim([min(min(XX)),max(max(XX))]);
ylim([min(min(XX)),max(max(XX))]);
caxis manual
caxis([min(min(twoDPot2DEG)) max(max(twoDPot2DEG))])
% color = [0.5,0.5,0.5];
color = 1;
ax = gca;
ax.Visible = 'off';

% Middle gate
gateMinZ = 17;
gateXWidth = 40;
gateYWidth = 40;
gateHoleHeight = 75;
metalThickness = 10;
gateOverHang = 7;

% for jj = 310
Xoffsets = [-120,-60,0,60,120];
% lightangle(h1,jj,50);
% lightangle(h2,90,-50);
for ii = 1:length(Xoffsets)
    
    % Hole part of the gate
    pts = [-gateXWidth/2 gateYWidth/2 gateMinZ;...
        -gateXWidth/2 -gateYWidth/2 gateMinZ;...
        gateXWidth/2 -gateYWidth/2 gateMinZ;...
        gateXWidth/2 gateYWidth/2 gateMinZ];
    fill3(pts(:,1) + Xoffsets(ii),pts(:,2),pts(:,3),color);
    pts = [-gateXWidth/2 gateYWidth/2 gateMinZ;...
        -gateXWidth/2 gateYWidth/2 gateMinZ + gateHoleHeight;...
        -gateXWidth/2 -gateYWidth/2 gateMinZ + gateHoleHeight;...
        -gateXWidth/2 -gateYWidth/2 gateMinZ];
    fill3(pts(:,1) + Xoffsets(ii),pts(:,2),pts(:,3),color);
    pts = [-gateXWidth/2 gateYWidth/2 gateMinZ;...
        -gateXWidth/2 gateYWidth/2 gateMinZ + gateHoleHeight;...
        gateXWidth/2 gateYWidth/2 gateMinZ + gateHoleHeight;...
        gateXWidth/2 gateYWidth/2 gateMinZ];
    fill3(pts(:,1) + Xoffsets(ii),pts(:,2),pts(:,3),color);
    pts = [gateXWidth/2 gateYWidth/2 gateMinZ;...
        gateXWidth/2 gateYWidth/2 gateMinZ + gateHoleHeight;...
        gateXWidth/2 -gateYWidth/2 gateMinZ + gateHoleHeight;...
        gateXWidth/2 -gateYWidth/2 gateMinZ];
    fill3(pts(:,1) + Xoffsets(ii),pts(:,2),pts(:,3),color);
    pts = [gateXWidth/2 -gateYWidth/2 gateMinZ;...
        gateXWidth/2 -gateYWidth/2 gateMinZ + gateHoleHeight;...
        -gateXWidth/2 -gateYWidth/2 gateMinZ + gateHoleHeight;...
        -gateXWidth/2 -gateYWidth/2 gateMinZ];
    fill3(pts(:,1) + Xoffsets(ii),pts(:,2),pts(:,3),color);

    % Arm part of the gate
    pts = [-gateXWidth/2 - gateOverHang max(max(YY)) gateMinZ + gateHoleHeight;...
        -gateXWidth/2 - gateOverHang -gateYWidth/2 - gateOverHang gateMinZ + gateHoleHeight;...
        gateXWidth/2 + gateOverHang -gateYWidth/2 - gateOverHang gateMinZ + gateHoleHeight;...
        gateXWidth/2 + gateOverHang max(max(YY)) gateMinZ + gateHoleHeight];
    fill3(pts(:,1) + Xoffsets(ii),pts(:,2),pts(:,3),color);
    pts = [-gateXWidth/2 - gateOverHang max(max(YY)) gateMinZ + gateHoleHeight + metalThickness;...
        -gateXWidth/2 - gateOverHang -gateYWidth/2 - gateOverHang gateMinZ + gateHoleHeight + metalThickness;...
        gateXWidth/2 + gateOverHang -gateYWidth/2 - gateOverHang gateMinZ + gateHoleHeight + metalThickness;...
        gateXWidth/2 + gateOverHang max(max(YY)) gateMinZ + gateHoleHeight + metalThickness];
    fill3(pts(:,1) + Xoffsets(ii),pts(:,2),pts(:,3),color);
    pts = [-gateXWidth/2 - gateOverHang max(max(YY)) gateMinZ + gateHoleHeight;...
        -gateXWidth/2 - gateOverHang max(max(YY)) gateMinZ + gateHoleHeight + metalThickness;...
        gateXWidth/2 + gateOverHang max(max(YY)) gateMinZ + gateHoleHeight +  metalThickness;...
        gateXWidth/2 + gateOverHang max(max(YY)) gateMinZ + gateHoleHeight];
    fill3(pts(:,1) + Xoffsets(ii),pts(:,2),pts(:,3),color);
    pts = [gateXWidth/2 + gateOverHang max(max(YY)) gateMinZ + gateHoleHeight;...
        gateXWidth/2 + gateOverHang max(max(YY)) gateMinZ + gateHoleHeight + metalThickness;...
        gateXWidth/2 + gateOverHang -gateYWidth/2 - gateOverHang gateMinZ + gateHoleHeight + metalThickness;...
        gateXWidth/2 + gateOverHang -gateYWidth/2 - gateOverHang gateMinZ + gateHoleHeight];
    fill3(pts(:,1) + Xoffsets(ii),pts(:,2),pts(:,3),color);
    pts = [gateXWidth/2 + gateOverHang -gateYWidth/2 - gateOverHang gateMinZ + gateHoleHeight;...
        gateXWidth/2 + gateOverHang -gateYWidth/2 - gateOverHang gateMinZ + gateHoleHeight + metalThickness;...
        -gateXWidth/2 - gateOverHang -gateYWidth/2 - gateOverHang gateMinZ + gateHoleHeight + metalThickness;...
        -gateXWidth/2 - gateOverHang -gateYWidth/2 - gateOverHang gateMinZ + gateHoleHeight];
    fill3(pts(:,1) + Xoffsets(ii),pts(:,2),pts(:,3),color);
    pts = [-gateXWidth/2 - gateOverHang -gateYWidth/2 - gateOverHang gateMinZ + gateHoleHeight;...
        -gateXWidth/2 - gateOverHang -gateYWidth/2 - gateOverHang gateMinZ + gateHoleHeight + metalThickness;...
        -gateXWidth/2 - gateOverHang max(max(YY)) gateMinZ + gateHoleHeight + metalThickness;...
        -gateXWidth/2 - gateOverHang max(max(YY)) gateMinZ + gateHoleHeight];
    fill3(pts(:,1) + Xoffsets(ii),pts(:,2),pts(:,3),color);
end
% drawnow;
% pause(0.05);
% end


%%
% Voltage pulse figure 5 gate
gT = [0,0.015,0.125,0.235,0.25,0.265,0.375,0.485,0.5,0.515,0.625,0.735,0.75,0.765,0.875,0.985,1.0]*8;
g1V = [0.8,0.8,0.8,0.791,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6];
g2V = [0.6,0.791,0.799,0.799,0.799,0.799,0.799,0.791,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6];
g3V = [0.6,0.6,0.6,0.6,0.6,0.791,0.7985,0.7985,0.7985,0.7985,0.7985,0.791,0.6,0.6,0.6,0.6,0.6];
g4V = [0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.791,0.799,0.799,0.799,0.799,0.799,0.791,0.6];
g5V = [0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.791,0.8,0.8,0.8];
figure;
hold on;
set(gca,'Fontsize',14);
plot(gT,g1V,'Linewidth',2.5);
plot(gT,g2V,'Linewidth',2.5);
plot(gT,g3V,'Linewidth',2.5);
plot(gT,g4V,'Linewidth',2.5);
plot(gT,g5V,'Linewidth',2.5);
xlabel('Time [ns]','Interpreter','Latex','Fontsize',22);
ylabel('Gate Voltage [V]','Interpreter','Latex','Fontsize',22);
% xlim([1,length(sparams.voltagePulse(1,:))]);
ylim([min(min(sparams.voltagePulse)),max(max(sparams.voltagePulse))*1.02]);
legend('V_1','V_2','V_3','V_4','V_5');


