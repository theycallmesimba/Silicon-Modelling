function checkLoeLOHOStates(sparams, XX, YY)
%CHECKLOELOHOSTATES Summary of this function goes here
%   Detailed explanation goes here
    for ii = 1:sparams.nSingleOrbitals
        for jj = 1:sparams.nSingleOrbitals
            currNorm = getInnerProduct2D(sparams.loeLOHOs(ii).wavefunctionMG,...
                sparams.loeLOHOs(jj).wavefunctionMG,XX,YY);
            % We will only print out the norm if it is below our acceptable
            % norm tolerance threshold.  Since we are doing things
            % numerically there is bound to be some error in our
            % orthogonality not being explicitly 0
            if ii == jj && abs(currNorm - 1) >= sparams.normThreshold
                fprintf(1,'Norm of Loe states not below threshold <%d|%d> = %g\n',ii,jj,currNorm);
            elseif ii ~= jj && abs(currNorm - 0) >= sparams.normThreshold
                fprintf(1,'Norm of Loe states not below threshold <%d|%d> = %g\n',ii,jj,currNorm);
            end   
        end
        if mod(ii-1,10) == 0
            fig = figure;
            pos = get(fig,'position');
            set(fig,'position',[pos(1:2)/4 pos(3)*1.8 pos(4)*1.3]);
        end
       
        subplot(2,5,mod(ii-1,10)+1);
        s = surf(XX,YY,sparams.loeLOHOs(ii).wavefunctionMG);
        set(s,'edgecolor','none');
        title(sprintf('%Leo LOHOs: %d',ii));
        view(2);
        colormap(hot);
        drawnow;
    end
end

