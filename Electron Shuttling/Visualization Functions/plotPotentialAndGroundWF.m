function fig = plotPotentialAndGroundWF( sparams, gateVoltages, xx )
%PLOTPOTENTIALANDGROUNDWF Summary of this function goes here
%   Detailed explanation goes here

    fig = figure;
    tempPot = sparams.P2DEGInterpolant([num2cell(gateVoltages),mat2cell(xx,1,length(xx))]);
    tempPot = squeezeFast(sparams.numOfGates,tempPot);
    
    [tempWF, ~] = solve1DSingleElectronSE(sparams,1,xx,tempPot); 
    xlabel('Position [nm]','Interpreter','Latex','Fontsize',14);
    xlim([min(xx),max(xx)]);
        
    yyaxis left
    plot(xx,tempPot/sparams.ee,'Linewidth',1.5);
    ylabel('Potential [V]','Interpreter','Latex','Fontsize',14);
    
    yyaxis right
    plot(xx,abs(tempWF).^2/norm(abs(tempWF).^2),'Linewidth',1.5);
    ylabel('Probability [Arb]','Interpreter','Latex','Fontsize',14);
    
    legend('Current Potential','Ground State |\psi|^2');
end

