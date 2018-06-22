% Test Ben's files....
v= [0.6,0.693,0.785,0.8];

for ii = 1:4
    for jj = 1:4
        for kk = 1:4
            for ll = 1:4
                for nn = 1:4
                    currFName = sprintf('V1_%.3f_V2_%.3f_V3_%.3f_V4_%.3f_V5_%.3f.csv',v(ii),v(jj),v(kk),v(ll),v(nn));
                    newFName = sprintf('NEWV1_%.3f_V2_%.3f_V3_%.3f_V4_%.3f_V5_%.3f.csv',v(kk),v(jj),v(ii),v(nn),v(ll));
                    movefile currFName newFName;
                end
            end
        end
    end
end
%%
currFName = 'V1_0.693_V2_0.600_V3_0.600_V4_0.800_V5_0.600.csv';

data = dlmread([sparams.potDir currFName]);
[rows,cols] = size(data);
xdata = data(1,2:cols);
ydata = data(2:rows,1);
zdata = data(2:rows,2:cols);

desiredGridX = 2^(nextpow2(length(xdata)));
desiredGridZ = round(2.5*length(ydata));

xxq = linspace(min(xdata),max(xdata),desiredGridX);
zzq = linspace(min(ydata),max(ydata),desiredGridZ);

[XX,ZZ] = meshgrid(xdata,ydata);
[XXq,ZZq] = meshgrid(xxq,zzq);

pot2D = -sparams.ee*interp2(XX,ZZ,zdata,XXq,ZZq); % Convert to J and invert
[~,sparams.twoDEGindZ] = min(abs(zz - (-0.5*1E-9)));
pot2DEG = pot2D(sparams.twoDEGindZ,:);

plot(xx*1E-9,pot2DEG)
% sparams.potentials(nn).gateValues = voltages(currFileVec);