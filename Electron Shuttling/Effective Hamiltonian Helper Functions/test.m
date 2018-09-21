profile on;
[~,nPts] = size(sparams.voltagePulse);

tc = zeros(1,nPts);
epsL = zeros(1,nPts);
epsR = zeros(1,nPts);

% for ii = [350,375,401,450]
for ii = 1:nPts
    if mod(ii,20) == 0
        fprintf(1,'(%04d/%04d)\n',ii,nPts);
    end
    currPotential = sparams.P2DEGInterpolant(getInterpolantArgument(sparams.voltagePulse(:,ii),xx));
    currPotential = squeezeFast(sparams.numOfGates,currPotential)';
    
    [tc(ii), epsL(ii), epsR(ii)] = calculateTunnelCoupling( sparams, xx, currPotential );
end
profile off
profile viewer
%%
epsL(epsL == 0) = NaN;
epsR(epsR == 0) = NaN;

figure;
yyaxis left
hold on;
plot(1:nPts,epsL/sparams.ee);
plot(1:nPts,epsR/sparams.ee);

yyaxis right
plot(sparams.voltagePulse');

figure;
tc(tc == 0) = NaN;
yyaxis left
hold on;
plot(1:nPts,tc/sparams.ee);

yyaxis right
plot(sparams.voltagePulse');
%%
% Query a particular point
index = [26];
vvs = zeros(length(index),length(xx));
for ii = 1:length(index)
    currPotential = sparams.P2DEGInterpolant(getInterpolantArgument(sparams.voltagePulse(:,index(ii)),xx));
    currPotential = squeezeFast(sparams.numOfGates,currPotential)';
    vvs(ii,:) = currPotential;
end

calculateTunnelCoupling(sparams, xx, vvs)/sparams.ee
figure;
hold on;
for ii = 1:length(index)
    plot(vvs(ii,:)/sparams.ee);
end

%%
aa = 30E-9;
omega = 1E13;
xxTemp = linspace(0,60E-9,200);
vvTemp = 1/2*sparams.me*omega^2*(xxTemp - aa).^2;
[ ~, ens ] = solve1DSingleElectronSE( sparams, 10, xxTemp, vvTemp );

ens/sparams.ee

eps0 = sparams.hbar*omega/2/sparams.ee
eps1 = sparams.hbar*omega/2*3/sparams.ee

%%
figure;
hold on;
plot(epsL);
plot(epsR);
legend('epsL','epsR');

figure;
plot(sparams.voltagePulse')

