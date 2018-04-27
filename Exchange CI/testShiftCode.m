xDom = linspace(-80,80,5000)*1E-9;
omega1 = 1.55E13;
omega2 = omega1;
m = 9.11E-31*0.191;
hbar = 6.626E-34/(2*pi);
a = 20E-9;

close all;
xx = X(1,:);
yy = Y(:,1);
ii = 3;
jj = 6;
omega = abs(mean(sparams.fittedPotentialParameters(:,1)));
alpha = sparams.me*omega/sparams.hbar;
xWF = 1/sqrt(2^ii*factorial(ii))*(alpha/pi)^(1/4)...
    *exp(-alpha*xx.^2/2).*hermiteH(ii,sqrt(alpha)*xx);
yWF = 1/sqrt(2^jj*factorial(jj))*(alpha/pi)^(1/4)...
    *exp(-alpha*yy.^2/2).*hermiteH(jj,sqrt(alpha)*yy);
figure;
plot(xx,xWF);
figure;
plot(yy,yWF);
figure;
WFMG = (yWF*ones(1,sparams.ngridx)).*(ones(sparams.ngridy,1)*xWF);
surf(X,Y,WFMG);
%%
% m = 1;
% hbar = 1;
a = 0;
% omega1 = 1;
% xDom = xDom/1E-9;

alpha1 = m*omega1/hbar;
alpha2 = m*omega2/hbar;
% Build shifted and non-shift potential wells
vorigin = 1/2*m*omega1^2*xDom.^2;
vshifted = 1/2*m*omega1^2*(xDom - a).^2;

% figure;
% hold on;
% plot(xDom, vorigin);
% plot(xDom, vshifted);
% for ii = 0:5
%     WF = 1/sqrt(2^ii*factorial(ii))*(alpha/pi)^(1/4)...
%         *exp(-alpha*xDom.^2/2).*hermiteH(ii,sqrt(alpha)*xDom);
%     WF = WF/norm(WF);
%     plot(xDom,WF);
% end

% Now, let's break down the shifted 6th HO in terms of unshifted HOs.
N = 6;
WForig = 1/sqrt(2^N*factorial(N))*(alpha1/pi)^(1/4)...
        *exp(-alpha1*(xDom - a).^2/2).*hermiteH(N,sqrt(alpha1)*(xDom - a));
norm(WForig)
WForig = WForig/norm(WForig);
norm(WForig)
plot(xDom,WForig)
%%
Nfit = 40;
bii = zeros(1,Nfit + 1);
for ii = linspace(0,Nfit,Nfit+1)
    fprintf(1,'Currently doing %d/%d\n',ii,Nfit);
%     tic;
    WFnonshifted = 1/sqrt(2^ii*factorial(ii))*(alpha2/pi)^(1/4)...
        *exp(-alpha2*xDom.^2/2).*hermiteH(ii,sqrt(alpha2)*xDom);
    WFnonshifted = WFnonshifted/norm(WFnonshifted);
    bii(ii+1) = sum(conj(WFnonshifted).*WForig);
%     toc;
end


% Check the decomposition
WFshifted = zeros(1,length(xDom));
for ii = linspace(0,Nfit,Nfit+1)
    WFnonshifted = 1/sqrt(2^ii*factorial(ii))*(alpha2/pi)^(1/4)...
        *exp(-alpha2*xDom.^2/2).*hermiteH(ii,sqrt(alpha2)*xDom);
    WFnonshifted = WFnonshifted/norm(WFnonshifted);
    
    WFshifted = WFnonshifted*bii(ii+1) + WFshifted;
    if mod(ii,2) == 0
        figure;
        plot(xDom,WFshifted);
        title(sprintf('Up to HO %d',ii));
    end
    
    % Calculate chi^2
    chi2 = sum((WFshifted - WForig).^2./WForig)/length(xDom);
    fprintf(1,'Norm up to %d is: %f\n', ii, norm(WFshifted));
    fprintf(1,'Current chi2 up to %d is: %f\n', ii, chi2);
end

figure;
plot(xDom,WForig);
title('Original shifted HO')




