function checkLOHOStates(sparams,XX,YY)
%CHECKLOHOSTATES Summary of this function goes here
%   Detailed explanation goes here
    for ii = 1:sparams.nSingleOrbitals
        currWFMG = sparams.LOHOs(ii).wavefunctionMG;
       
        if mod(ii-1,8) == 0
            fig = figure;
            pos = get(fig,'position');
            set(fig,'position',[pos(1:2)/4 pos(3)*1.8 pos(4)*1.3]);
        end
        subplot(2,4,mod(ii-1,8)+1);
        
        s = surf(XX,YY,currWFMG);
        set(s,'edgecolor','none');
        title(sprintf('Dot: %d  N: %d  M: %d',sparams.LOHOs(ii).dot,...
            sparams.LOHOs(ii).n,sparams.LOHOs(ii).m));
        view(2);
        colormap(jet);
       
        fprintf(1,'LOHO state %03d - dot: %d, n: %d, m: %d, energy: %0.3E [J], norm: %0.3E \n',...
            ii, sparams.LOHOs(ii).dot, sparams.LOHOs(ii).n,sparams.LOHOs(ii).m,...
            sparams.LOHOs(ii).energy, getInnerProduct2D(currWFMG,currWFMG,XX,YY));
        drawnow;
    end
end

