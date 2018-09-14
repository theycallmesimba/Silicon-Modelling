function checkOriginHOStates(sparams, XX, YY)
%CHECKORIGINHOSTATES Summary of this function goes here
%   Detailed explanation goes here
    for ii = 1:sparams.nOriginHOs
        currWFMG = sparams.originHOs(ii).wavefunctionMG;
       
        if mod(ii-1,8) == 0
            fig = figure;
            pos = get(fig,'position');
            set(fig,'position',[pos(1:2)/4 pos(3)*1.8 pos(4)*1.3]);
        end
        subplot(2,4,mod(ii-1,8)+1);
        
        s = surf(XX,YY,currWFMG);
        set(s,'edgecolor','none');
        title(sprintf('Origin HOs N: %d  M: %d',...
            sparams.originHOs(ii).n,sparams.originHOs(ii).m));
        view(2);
        colormap(jet);
       
        fprintf(1,'Origin HO state %03d - n: %d, m: %d, energy: %0.3E [J], norm: %0.3E \n',...
            ii, sparams.originHOs(ii).n,sparams.originHOs(ii).m,...
            sparams.originHOs(ii).energy, getInnerProduct2D(currWFMG,currWFMG,XX,YY));
        drawnow;
    end
end

