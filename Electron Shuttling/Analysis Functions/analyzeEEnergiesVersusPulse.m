function analyzeEEnergiesVersusPulse( sparams, xx )
%ANALYZEEIEGENENERGIESVERSUSPULSE Summary of this function goes here
%   Detailed explanation goes here

    sparams.vPulseGInterpolants = {};
    xPots = linspace(0,1,length(sparams.voltagePulse(1,:)));    

    for vv = 1:sparams.numOfGates
        vPulseGInterpolants{vv} = griddedInterpolant({xPots},sparams.voltagePulse(vv,:));
    end
    
    qPots = linspace(0,1,301);
    
    pulse = getInterpolatedPulseValues(sparams,qPots,vPulseGInterpolants);
    
    nSols = 4;
    ens = zeros(length(qPots),nSols);
    for ii = 1:length(qPots)
        currPot = sparams.P2DEGInterpolant(getInterpolantArgument(pulse(:,ii),xx));
        [~,energies] = solve1DSingleElectronSE(sparams, nSols, xx, currPot);
        
        ens(ii,:) = diag(energies)'/sparams.ee;
    end
    
    figure;
    axisFontSize = 20;
    tickFontSize = 15;
    titleFontSize = 25;
    set(gca,'Fontsize',tickFontSize);
    hold on;
    title('E-Energies vs. Pulse','Interpreter','Latex','Fontsize',titleFontSize);
    xlabel('\% of shuttle time','Interpreter','Latex','Fontsize',axisFontSize);
    
    yyaxis left
    plot(qPots,ens,'Linewidth',2.0,'Linestyle','-');
    ylabel('Energy [eV]','Interpreter','Latex','Fontsize',axisFontSize);
    
    yyaxis right
    for ii = 1:sparams.numOfGates
        plt = plot(xPots,sparams.voltagePulse(ii,:),'Linewidth',2.0,...
            'Linestyle','-','Marker','none');
        plt.Color(4) = 0.7;
    end
    ylabel('Voltage [V]','Interpreter','Latex','Fontsize',axisFontSize);
end

