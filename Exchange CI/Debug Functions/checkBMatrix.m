function checkBMatrix(sparams,XX,YY)
%CHECKBMATRIX Summary of this function goes here
%   Detailed explanation goes here

    for jj = 1:sparams.nSingleOrbitals
        currWF = zeros(sparams.ngridy,sparams.ngridx);
        for kk = 1:sparams.nOriginHOs
            currWF = currWF + sparams.bcoeffs(jj,kk)*sparams.originHOs(kk).wavefunctionMG;
        end

        fig = figure;
        pos = get(fig,'position');
        set(fig,'position',[pos(1:2)/4 pos(3)*1.8 pos(4)*1.3]);
        subplot(1,2,1);
        s = surf(XX,YY,currWF);
        title(sprintf('B transformed %d',jj));
        set(s,'edgecolor','none');
        colormap(hot)
        view(2);
        
        subplot(1,2,2);
        s = surf(XX,YY,sparams.loeLOHOs(jj).wavefunctionMG);
        set(s,'edgecolor','none');
        title(sprintf('LoeHO %d',jj));
        view(2);
        colormap(hot);
        drawnow;

        fprintf(1,'State %03d norms - b transformed: %0.4f, LoeHO states: %0.4f \n',...
            jj, getInnerProduct2D(currWF,currWF,XX,YY),...
            getInnerProduct2D(sparams.loeLOHOs(jj).wavefunctionMG,...
            sparams.loeLOHOs(jj).wavefunctionMG,XX,YY));
    end
end

