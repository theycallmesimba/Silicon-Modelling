function plotEffHamResults2D(dataX, dataStrX, dataY, dataStrY, dataZ, dataStrZ)
%PLOTEFFHAMRESULTS Summary of this function goes here
%   Detailed explanation goes here

    axisFontSize = 15;
    labelFontSize = 20;
    
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
        case 'dPhi'
            xLabelStr = 'Valley splitting phase difference [rads]';
            xInd = 5;
        case 'spinOrbit'
            xLabelStr = 'Spin orbit coupling [$\mu$eV]';
            dataX = dataX/1.602E-19/1E-6;
            xInd = 6;
        case 'BfieldZ'
            xLabelStr = 'Magnetic field $z$-component [meV]';
            xInd = 7;
        case 'BfieldX'
            xLabelStr = 'Magnetic field $x$-component [$\mu$eV]';
            xInd = 8;
        otherwise
            error('Could not find matching string for X axis variable.');
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
        case 'dPhi'
            yLabelStr = 'Valley splitting phase difference [rads]';
            yInd = 5;
        case 'spinOrbit'
            yLabelStr = 'Spin orbit coupling [$\mu$eV]';
            dataY = dataY/1.602E-19/1E-6;
            yInd = 6;
        case 'BfieldZ'
            yLabelStr = 'Magnetic field $z$-component [meV]';
            yInd = 7;
        case 'BfieldX'
            yLabelStr = 'Magnetic field $x$-component [$\mu$eV]';
            yInd = 8;
        otherwise
            error('Could not find matching string for Y axis variable.');
    end
    
    switch dataStrZ
        case 'velocity'
            zAxisLabel = 'Electron velocity [nm/ns]';
        case 'pulseTime'
            zAxisLabel = 'Pulse time [s]';
        case 'fidelity'
            zAxisLabel = 'Final state fidelity';
        otherwise
            error('Could not find matching string for Y axis variable.');
    end
    zAxisLabel
    
    % Just in case use forgot to squeeze before sending data...
    dataZ = squeeze(dataZ);
    
    if yInd < xInd
        DZ = dataZ;
    elseif yInd == xInd
        fprintf(1,'X and Y data are the same variable.\n');
        return;
    else
        DZ = dataZ';
    end
   
    figure;
    hold on;
    [DX,DY] = meshgrid(dataX,dataY);
    s = surf(DX,DY,DZ);
    
    if max(dataX)/min(dataX) >= 100
        set(gca,'XScale','log');
    end
    if max(dataY)/min(dataY) >= 100
        set(gca,'YScale','log');
    end
    
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

