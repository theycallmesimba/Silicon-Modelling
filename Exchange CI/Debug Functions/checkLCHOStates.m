function checkLCHOStates(sparams, XX, YY)
%CHECKLCHOSTATES Summary of this function goes here
%   Detailed explanation goes here

    for ii = 1:sparams.nSingleOrbitals     
        currWF = sparams.LCHOs(ii).wavefunctionMG;

        if mod(ii-1,8) == 0
            fig = figure;
            pos = get(fig,'position');
            set(fig,'position',[pos(1:2)/4 pos(3)*1.8 pos(4)*1.3]);
        end
        subplot(2,4,mod(ii-1,8)+1);
        s = surf(XX,YY,currWF);
        set(s,'edgecolor','none');
        title(sprintf('LCHO %d',ii));
        view(2);
        colormap(jet);

        fprintf(1,'LCHO state %03d - energy: %0.3E [J]\n',...
            ii, sparams.LCHOs(ii).energy);
    end
end

