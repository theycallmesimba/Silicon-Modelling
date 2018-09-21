function [purityZZ, T2Overtc] = chargeDecoherenceMain( sparams, vPulse, vPulseTime, xx )
%CHARGEDECOHERENCEMAIN Summary of this function goes here
%   Detailed explanation goes here
    
    T2Sweep = 100000E-9;
    
    pulseTime = vPulseTime;
    
    purityZZ = zeros(length(T2Sweep));
    T2Overtc = zeros(length(T2Sweep));
    
    [epsL, epsR] = getDetuningVsVoltagePulse(sparams, xx, vPulse);
    
    mm = 0;
    for jj = 1:length(T2Sweep)
        tic;

        mm = mm + 1;

        T2 = T2Sweep(jj);

        fprintf(1,'(%d/%d): T2 = %.3E [s]\n',mm,length(T2Sweep),T2);

        pulseTVec = linspace(0,pulseTime,length(vPulse(1,:)));
        pulseInterpolants = makePulseInterpolants(sparams,pulseTVec,vPulse);
        
        detInterpolantL = griddedInterpolant({pulseTVec},epsL);
        detInterpolantR = griddedInterpolant({pulseTVec},epsR);
        
        dt = 1E-15;
        time = 0:dt:pulseTime;
%         finePulse = getInterpolatedPulseValues(sparams,time,pulseInterpolants);
        fineEpsL = detInterpolantL({time});
        fineEpsR = detInterpolantR({time});
        
        tcVArg = vPulse(:,1);
        tcVArg(sparams.gatesUsedInPulse(2)) = max(vPulse(sparams.gatesUsedInPulse(2),:));
        currPotential = sparams.P2DEGInterpolant(getInterpolantArgument(tcVArg,xx));
        currPotential = squeezeFast(sparams.numOfGates,currPotential)';
        tc = calculateTunnelCoupling( sparams, xx, currPotential );
        
        h = waitbar(0,sprintf('Current Time Index: %d/%d',0,length(time)),...
            'Name','Simulating charge decoherence...',...
            'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');

%         rho0 = [1,0;0,0];
        H0 = [fineEpsL(1),tc;tc,fineEpsR(1)];
        [kets,ens] = eig(H0);
        [~, ind] = sort(diag(ens));
        ens = ens(ind,ind);
        kets = kets(:,ind);
%         kets(:,1)
        
        currRho = kets(:,1)*kets(:,1).';

        expectationValue = zeros(3,length(time));
        expectationValue(1,1) = trace(currRho*sparams.dz);
        expectationValue(2,1) = trace(currRho*sparams.dx);
        expectationValue(3,1) = trace(currRho*sparams.dy);
        purity = zeros(1,length(time));
        purity(1) = trace(currRho*currRho);
        for ii = 2:length(time)
            % Check for cancel button click
            if getappdata(h,'canceling')
                break;
            end

            % Update waitbar every N frames
            if mod(ii,4000) == 0
                waitbar(ii/length(time), h, sprintf('Current Time Index: %d/%d',ii,length(time)));
            end
            
%             currPotential = sparams.P2DEGInterpolant(getInterpolantArgument(finePulse(:,ii),xx));
%             currPotential = squeezeFast(sparams.numOfGates,currPotential)';
%             [tc, epsL, epsR] = calculateTunnelCoupling(sparams, xx, currPotential);
            
            H = [fineEpsL(ii), tc;tc, fineEpsR(ii)];
            currRho = currRho + dt*getCurrentMasterEquationOperator(sparams,currRho,H,T2);
            expectationValue(1,ii) = trace(currRho*sparams.dz);
            expectationValue(2,ii) = trace(currRho*sparams.dx);
            expectationValue(3,ii) = trace(currRho*sparams.dy);
            purity(ii) = trace(currRho*currRho);
        end
        [kets,ens] = eig(H);
        [~, ind] = sort(diag(ens));
        ens = ens(ind,ind);
        kets = kets(:,ind);
        finalRho0 = kets(:,1)*kets(:,1).';
        finalFid = trace(finalRho0*currRho)
        
        delete(h);
        purityZZ(jj) = purity(end);
        toc;
    end
    
    figure;
    plot(time,expectationValue,'Linewidth',2.0);
    ylim([-1.1,1.1]);
    ylabel('$\langle A\rangle$','Interpreter','latex','Fontsize',20);
    xlabel('Time [s]','Interpreter','latex','Fontsize',20);
    legend('<Z>','<X>','<Y>');
    
    figure;
    plot(time,purity,'Linewidth',2.0);
    ylim([0,1.1]);
    ylabel('${\rm Tr} [\rho^2]$','Interpreter','latex','Fontsize',20);
    xlabel('Time [s]','Interpreter','latex','Fontsize',20);
    
end

function dRho_dt = getCurrentMasterEquationOperator(sparams,rho,H,T2)
    unitaryTerm = -1i/sparams.hbar*(H*rho - rho*H);
    
    linOp = sqrt(1/T2)*[0,0;0,1];
    dissipatorTerm = -1/2*((linOp')*linOp*rho + rho*(linOp')*linOp) + linOp*rho*linOp;
    
    dRho_dt = unitaryTerm + dissipatorTerm;
end













