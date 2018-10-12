% Zeroth step, load the potentials
shuttleParameterFile;

fprintf(1,'Loading potentials...\n');
[sparams,xx,zz] = loadPotentials(sparams);
sparams.nxGrid = length(xx);
sparams.nzGrid = length(zz);
sparams.dx = xx(2) - xx(1);
sparams.dz = zz(2) - zz(1);

% Find which index corresponds to where the 2DEG should be
[~,sparams.twoDEGindZ] = min(abs(zz - (-0.5*1E-9)));
for ii = 1:length(sparams.potentials)
    sparams.potentials(ii).pot2DEG = sparams.potentials(ii).pot2D(sparams.twoDEGindZ,:);
end

% Now we want to make the potential interpolant object (both 2D and 2DEG)
sparams.interpType = 'spline';
sparams.extrapType = 'spline';
sparams = makePotentialsInterpolants(sparams,xx,zz);

% First step is to define our voltage search range and then find how the
% detuning changes in that window.
adiabThresh = [0.01,0.01];
vBounds = [0.44,0.60; 0.44,0.6; 0.4,0.6];

[~, voltagePulse, ~] = getVoltagePulseAdiabatic( sparams, xx, adiabThresh, vBounds, 0, NaN );
[epsL, epsR] = getDetuningVsVoltagePulse( sparams, xx, voltagePulse, 1 );

% figure;
% hold on;
% plot(epsL);
% plot(epsR);
% 
% figure;
% hold on;
% for ii = 1:sparams.numOfGates
%     plot(voltagePulse(ii,:));
% end

% Second step is to get the tunnel coupling at the maximum voltage point in
% our pulse
tcVArg = voltagePulse(:,1);
tcVArg(sparams.gatesUsedInPulse(2)) = max(voltagePulse(sparams.gatesUsedInPulse(2),:));
currPotential = sparams.P2DEGInterpolant(getInterpolantArgument(tcVArg,xx));
currPotential = squeezeFast(sparams.numOfGates,currPotential)';
tc = calculateTunnelCoupling( sparams, xx, currPotential );

tc/sparams.ee
%%
% Third step is to sweep over spin orbit, valley, T2, etc and calculate an
% adiabatic pulse
% profile on
T2Sweep = 100000E-9;
adiabThresh = [0.005];
vBounds = [0.44,0.60; 0.44,0.6; 0.4,0.6];
sparams.dt = 5E-14;

dPhi = [0]*pi;
valleyL = [30,60,90]*1E-6*sparams.ee;
valleyR = [30,60,90]*1E-6*sparams.ee;
spinOrbit = [0.2,1.0,2.0]*1E-6*sparams.ee;
Ez = [40]*1E-6*sparams.ee;
Ex = [1.6]*1E-6*sparams.ee;
% dphi = 0.132*pi;
% valleyL = [0]*1E-6*exp(1i*dphi)*sparams.ee;
% valleyR = [0]*1E-6*sparams.ee;
% spinOrbit = [0]*1E-6*sparams.ee;
% Ez = [0]*1E-6*sparams.ee;
% Ex = [0]*1E-6*sparams.ee;

pTime = zeros(length(adiabThresh),length(T2Sweep),length(valleyL),length(valleyR),...
    length(dPhi),length(spinOrbit),length(Ez),length(Ex));
fidelity = zeros(length(adiabThresh),length(T2Sweep),length(valleyL),length(valleyR),...
    length(dPhi),length(spinOrbit),length(Ez),length(Ex),sparams.nFidelityFrames);
purity = zeros(length(adiabThresh),length(T2Sweep),length(valleyL),length(valleyR),...
    length(dPhi),length(spinOrbit),length(Ez),length(Ex),sparams.nPurityFrames);
orbExp = zeros(length(adiabThresh),length(T2Sweep),length(valleyL),length(valleyR),...
    length(dPhi),length(spinOrbit),length(Ez),length(Ex),3,sparams.nExpectationFrames);
orbFRho = zeros(length(adiabThresh),length(T2Sweep),length(valleyL),length(valleyR),...
    length(dPhi),length(spinOrbit),length(Ez),length(Ex),2,2);
valExp = zeros(length(adiabThresh),length(T2Sweep),length(valleyL),length(valleyR),...
    length(dPhi),length(spinOrbit),length(Ez),length(Ex),3,sparams.nExpectationFrames);
valFRho = zeros(length(adiabThresh),length(T2Sweep),length(valleyL),length(valleyR),...
    length(dPhi),length(spinOrbit),length(Ez),length(Ex),2,2);
spinExp = zeros(length(adiabThresh),length(T2Sweep),length(valleyL),length(valleyR),...
    length(dPhi),length(spinOrbit),length(Ez),length(Ex),3,sparams.nExpectationFrames);
spinFRho = zeros(length(adiabThresh),length(T2Sweep),length(valleyL),length(valleyR),...
    length(dPhi),length(spinOrbit),length(Ez),length(Ex),2,2);

uu = 1;
mm = 0;
for cc = 1:length(adiabThresh)
    for ii = 1:length(T2Sweep)
        for jj = 1:length(valleyL)
        for kk = 1:length(valleyR)
        for pp = 1:length(dPhi) 
            % If we don't need valley, just omit from the simulation as it
            % messes with the adiabatic pulse finding
            if valleyL(jj) == 0 || valleyR(kk) == 0
                sparams.includeValley = 0;
            else
                sparams.includeValley = 1;
            end

            for aa = 1:length(spinOrbit)
                for bb = 1:length(Ez)
                for dd = 1:length(Ex)
                    % If we don't need spin, just omit it from the
                    % simulation as it messes with the adiabatic pulse
                    % finding
                    if spinOrbit(aa) == 0 || Ez(bb) == 0 || Ex(dd) == 0
                        sparams.includeSpin = 0;
                    else
                        sparams.includeSpin = 1;
                    end

                    % Display to use what parameters are bein probed right
                    % now
                    mm = mm + 1;
                    fprintf(1,'(%d/%d): Adiabatic Threshold = %.3E [arb]\n',...
                        mm,length(T2Sweep)*length(valleyL)*length(valleyR)*length(dPhi)*...
                        length(spinOrbit)*length(Ez)*length(Ex)*length(spinOrbit)*length(adiabThresh),...
                        adiabThresh(cc));
                    fprintf(1,'T2 = %.3E [s]\n',T2Sweep(ii));
                    fprintf(1,'ValleyL = %.3E [eV], ValleyR = %.3E [eV], Phase = %.3f\n',...
                        valleyL(jj)/sparams.ee, valleyR(kk)/sparams.ee, dPhi(pp));
                    fprintf(1,'Ez = %.3E [eV], Ex = %.3E [eV], Spin Orbit = %.3E [eV]\n',...
                        Ez(bb)/sparams.ee, Ex(dd)/sparams.ee, spinOrbit(aa)/sparams.ee);

                    % Put all the hamiltonian parameters into a single
                    % variable for ease
                    effHamiltonianParams = buildEffHamiltonianParamVariable(epsL, epsR,...
                        tc, Ez(bb), Ex(dd), valleyL(jj)*exp(1i*dPhi(pp)), valleyR(kk), spinOrbit(aa), spinOrbit(aa));

                    % Get an adiabatic voltage pulse based on the voltage
                    % points given earlier and using the effective
                    % Hamiltonian parameters
                    [sparams, voltagePulse, pTimeTemp] = getVoltagePulseAdiabatic(...
                        sparams, xx, [adiabThresh(cc),adiabThresh(cc)], vBounds, 1, effHamiltonianParams );
                    fprintf(1,'Pulse found!\nPulse time = %.6E [s]\n',pTimeTemp);

                    % Using the adiabatic voltage pulse, run an effective
                    % shuttling simulation.
                    fprintf(1,'Now doing effective shuttling simulation using found pulse...\n\n');
                    [fidTemp, purTemp, orbExpTemp, valExpTemp, spinExpTemp, fRho] =...
                        simulateEffectiveShuttling(sparams, xx, voltagePulse, pTimeTemp, effHamiltonianParams, T2Sweep(ii));

                    % Store some data and calculate some partial traces
                    pTime(cc,ii,jj,kk,pp,aa,bb,dd) = pTimeTemp;
                    fidelity(cc,ii,jj,kk,pp,aa,bb,dd,:) = fidTemp;
                    purity(cc,ii,jj,kk,aa,bb,dd,:) = purTemp;
                    orbExp(cc,ii,jj,kk,pp,aa,bb,dd,:,:) = orbExpTemp;
                    orbFRho(cc,ii,jj,kk,pp,aa,bb,dd,:,:) = partialTrace(fRho,[2,3],[2,2,2]); 
                    valExp(cc,ii,jj,kk,pp,aa,bb,dd,:,:) = valExpTemp;
                    valFRho(cc,ii,jj,kk,pp,aa,bb,dd,:,:) = partialTrace(fRho,[1,3],[2,2,2]); 
                    spinExp(cc,ii,jj,kk,pp,aa,bb,dd,:,:) = spinExpTemp;
                    spinFRho(cc,ii,jj,kk,pp,aa,bb,dd,:,:) = partialTrace(fRho,[1,2],[2,2,2]);   
                end
                end
            end
        end
        end 
        end
    end
end
% profile viewer
% profile off
%%
nPts = 500;
dPhi = 0.132*pi;
valleyL = [60]*1E-6*exp(1i*dPhi)*sparams.ee;
valleyR = [30]*1E-6*sparams.ee;
spinOrbit = [0.2]*1E-6*sparams.ee;
Ez = [40]*1E-6*sparams.ee;
Ex = [1.6]*1E-6*sparams.ee;
tc = 45E-6*sparams.ee;

effHamiltonianParams = buildEffHamiltonianParamVariable(NaN, NaN,...
    tc, Ez, Ex, valleyL, valleyR, spinOrbit, spinOrbit);

epsL = linspace(-0.2,0.2,nPts)*1E-3*sparams.ee;
epsR = linspace(0.2,-0.2,nPts)*1E-3*sparams.ee;

energies = zeros(8,nPts);
for ii = 1:nPts
    effHamiltonianParams{1} = epsL(ii);
    effHamiltonianParams{2} = epsR(ii);
    Heff = constructEffectiveHamiltonian( sparams, effHamiltonianParams);
    
    [kets,ens] = eig(Heff);
    [~, ind] = sort(diag(ens));
    ens = ens(ind,ind);
    energies(:,ii) = diag(ens);
end

plot(epsL,energies/sparams.ee,'Linewidth',2);
%%
sparams.includeValley = 1;
effHamParams = buildEffHamiltonianParamVariable(0, 0, tc, NaN, valleySplit(jj), valleySplit(jj), NaN, NaN);
effHamTemp = constructEffectiveHamiltonian( sparams, effHamParams);
[ketsph,ensph] = eig(effHamTemp);
ensph = ensph/sparams.ee;
[~, ind] = sort(diag(ensph));
ketsph = ketsph(:,ind);
%%
