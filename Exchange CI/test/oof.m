T2 = 3E-9;
t = logspace(-13,-9,300);

P = 1/2*(1 + exp(-2*t/T2));
F = 1/2*(1 + exp(-t/T2));

figure;
loglog(t,1-P,t,1-F);
grid on
grid minor

%%
error = linspace(0.6,0.99,300);
nSteps = 30;
T2 = 3E-9;

P = error.^(1/nSteps);
t = -log(2*P - 1)/2*T2;

plot(error,t);
set(gca,'Fontsize',14);
xlabel('Singlet distribution total error','Interpreter','latex','Fontsize',22);
ylabel('Time per shuttle','Interpreter','latex','Fontsize',22);
ylim([min(t),max(t)]);
line([0.9 0.9], [min(t) max(t)],'Color','r','Linewidth',2);
grid minor
grid on