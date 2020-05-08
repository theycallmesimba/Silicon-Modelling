function checkBasisTransformation( sparams, gridparams, basis1, basis2, Tmatrix )
%CHECKBASISTRANSFORMATION Summary of this function goes here
%   Detailed explanation goes here
    [~,size1] = size(basis1);
    [~,size2] = size(basis2);
    
    XX = gridparams.XX;
    YY = gridparams.YY;

    nStatesToCompare = min([size1,size2]);
    maxPlotsPerFigure = 8;
    foundNonNorm = 0;
    
    for ii = 1:nStatesToCompare
        tempwf = zeros(gridparams.ngridy,gridparams.ngridx);

        for jj = 1:size1
            tempwf = tempwf + Tmatrix(ii,jj)*basis1(jj).wavefunctionMG;
        end
        
        currNorm = getInnerProduct2D(tempwf,tempwf,XX,YY);
        
        if abs(currNorm - 1) >= sparams.normThreshold
            foundNonNorm = 1;
            fprintf(1,'WARNING: Transformed state |%d> is not normalized: %g\n',ii,currNorm);
        end
        
        % Now plot the transformed wave function along with the
        % wavefunction we wish to compare with from basis2.
        if mod(2*(ii-1),maxPlotsPerFigure) == 0
            fig = figure;
            pos = get(fig,'position');
            set(fig,'position',[pos(1:2)/4 pos(3)*1.8 pos(4)*1.3]);
        end
        subplot(2,maxPlotsPerFigure/2,mod(2*ii-1,maxPlotsPerFigure));
        s = surf(XX,YY,tempwf);
        title(sprintf('Basis 1 Transformed: %d',ii));
        set(s,'edgecolor','none');
        xlim([min(min(XX)),max(max(XX))]);
        ylim([min(min(YY)),max(max(YY))]);
        colormap(hot)
        view(2);

        pltNum = mod(2*ii,maxPlotsPerFigure);
        if pltNum == 0
            pltNum = maxPlotsPerFigure;
        end
        subplot(2,maxPlotsPerFigure/2,pltNum);
        s = surf(XX,YY,basis2(ii).wavefunctionMG);
        set(s,'edgecolor','none');
        title(sprintf('Basis 2: %d',ii));
        xlim([min(min(XX)),max(max(XX))]);
        ylim([min(min(YY)),max(max(YY))]);
        view(2);
        colormap(hot);
        drawnow;
    end
    
    if ~foundNonNorm
        fprintf(1,'All transformed states are normalized.\n');
    end
end

