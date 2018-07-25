function analyzeStarkShift( sparams )
%ANALYZESTARKSHIFT Summary of this function goes here
%   Detailed explanation goes here

    sweepVec = getSweepVector(sparams);
    
    if sparams.calculateStarkShift
        for ii = 1:length(sweepVec)
            if mod(ii-1,6) == 0
                fig = figure;
                pos = get(fig,'position');
                set(fig,'position',[pos(1:2)/4 pos(3)*2.0 pos(4)*1.25]);
            end

            subplot(2,3,mod(ii-1,6)+1);
            hold on;
            plot(sparams.tStarkShift(ii,:),(sparams.vShift(ii,:)-sparams.vShift(ii,1))*1E-6);
            plot(sparams.tStarkShift(ii,:),(sparams.vShiftGround(ii,:)-sparams.vShift(ii,1))*1E-6);
            xlabel('Time index [s]','interpreter','latex','fontsize',10);
            ylabel('$\nu - \nu_0$ [MHz]','interpreter','latex','fontsize',10);
            title(['Shuttling Simulation ' num2str(sparams.totalTime(ii)) '[s]'],'interpreter','latex','fontsize',10);
            drawnow;
        end
    end
end

