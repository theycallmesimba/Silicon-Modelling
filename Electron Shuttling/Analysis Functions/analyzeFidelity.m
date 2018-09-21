function analyzeFidelity( sparams )
%ANALYZEFIDELITY Summary of this function goes here
%   Detailed explanation goes here

    sweepVec = getSweepVector(sparams);
    fids = sparams.fidelity;
    
    figure;
    plot(linspace(0,1,length(fids(1,:))),fids','Linewidth',1.5);
    set(gca,'Fontsize',14);
    xlabel('t/T','Fontsize',22,'Interpreter','Latex');
    ylabel('Fidelity','Fontsize',22,'Interpreter','Latex');
    
    % Make waitbar
    if strcmp(sparams.sweptParameter,'time')        
        legVec = cell(1,length(sweepVec));
        for ii = 1:length(sweepVec)
            legVec{ii} = num2str(sweepVec(ii));
        end
    elseif strcmp(sparams.sweptParameter,'adiabicity')     
        [rows,~] = size(sweepVec);
        legVec = cell(1,rows);
        for ii = 1:rows
            legVec{ii} = [num2str(sweepVec(ii,1)) '-' num2str(sweepVec(ii,2))];
        end
    end
    legend(legVec);
end

