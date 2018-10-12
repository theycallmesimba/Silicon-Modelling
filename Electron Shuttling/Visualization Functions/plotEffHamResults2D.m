function plotEffHamResults2D(dataX, dataY, dataStrX, dataStrY, effHamResultsData)
%PLOTEFFHAMRESULTS Summary of this function goes here
%   Detailed explanation goes here

    axisFontSize = 15;
    labelFontSize = 20;
    
    zAxisLabel = 'Electron velocity [nm/ns]';
    
    % Just in case use forgot to squeeze before sending data...
    effHamResultsData = squeeze(effHamResultsData);
    
    switch dataStrX
        case 'adiabThresh'
            xLabelStr = 'Adiabatic parameter [arb]';
            xInd = 1;
        case 'T2Sweep'
            xLabelStr = 'T_2 [s]';
            xInd = 2;
        case 'valleyL'
            xLabelStr = 'Valley splitting (left dot) [$\mu$eV]';
            dataX = dataX/1.602E-19/1E-6;
            xInd = 3;
        case 'valleyR'
            xLabelStr = 'Valley splitting (right dot) [$\mu$eV]';
            dataX = dataX/1.602E-19/1E-6;
            xInd = 4;
        case 'spinOrbit'
            xLabelStr = 'Spin orbit coupling [$\mu$eV]';
            dataX = dataX/1.602E-19/1E-6;
            xInd = 5;
        case 'Bfield'
            xLabelStr = 'Magnetic field $z$-component [T]';
            xInd = 6;
        otherwise
            fprintf(1,'Could not find matching string for X axis variable.\n');
            return;
    end
    
    switch dataStrY
        case 'adiabThresh'
            yLabelStr = 'Adiabatic parameter [arb]';
            yInd = 1;
        case 'T2Sweep'
            yLabelStr = 'T_2 [s]';
            yInd = 2;
        case 'valleyL'
            yLabelStr = 'Valley splitting (left dot) [$\mu$eV]';
            dataY = dataY/1.602E-19/1E-6;
            yInd = 3;
        case 'valleyR'
            yLabelStr = 'Valley splitting (right dot) [$\mu$eV]';
            dataY = dataY/1.602E-19/1E-6;
            yInd = 4;
        case 'spinOrbit'
            yLabelStr = 'Spin orbit coupling [$\mu$eV]';
            dataY = dataY/1.602E-19/1E-6;
            yInd = 5;
        case 'Bfield'
            yLabelStr = 'Magnetic field $z$-component [T]';
            yInd = 6;
        otherwise
            fprintf(1,'Could not find matching string for Y axis variable.\n');
            return;
    end
    
    if yInd < xInd
        DZ = effHamResultsData;
    elseif yInd == xInd
        fprintf(1,'X and Y data are the same variable.\n');
        return;
    else
        DZ = effHamResultsData';
    end
   
    figure;
    hold on;
    [DX,DY] = meshgrid(dataX,dataY);
    s = surf(DX,DY,DZ);
    
    set(gca,'Fontsize',axisFontSize);
    set(s,'EdgeColor','none');
    cb = colorbar;
    cb.Label.String = zAxisLabel;
    cb.Label.FontSize = labelFontSize;
    cb.Label.Interpreter = 'latex';
    xlabel(xLabelStr,'Fontsize',labelFontSize,'Interpreter','latex');
    ylabel(yLabelStr,'Fontsize',labelFontSize,'Interpreter','latex');
    xlim([min(dataX), max(dataX)]);
    ylim([min(dataY), max(dataY)]);
end

