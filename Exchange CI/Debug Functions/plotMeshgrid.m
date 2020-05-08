function plotMeshgrid( gridparams, VV )
%PLOTMESHGRID Summary of this function goes here
%   Detailed explanation goes here

    figure('Color','white');
    surf(gridparams.XX,gridparams.YY,VV);
    set(gca,'Fontsize',14,'TicklabelInterpreter','latex');
    xlabel('x-axis','Fontsize',18,'Interpreter','latex');
    ylabel('y-axis','Fontsize',18,'Interpreter','latex');
    xlim([min(min(gridparams.XX)),max(max(gridparams.XX))]);
    ylim([min(min(gridparams.YY)),max(max(gridparams.YY))]);
    colormap(viridis());
    shading interp
    view(2);
end

