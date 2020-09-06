Initialize_Potential;

% Simulation parameters
% Dot potential size along x (characteristic width)
% lxi = 2.3;
% lxi = [1.5,2,2.5,3.0,3.5];
lxi = linspace(1.5,3.5,31);
% Potential minima
Vimin = 10;
% Dot separation (x-axis)
dotSep = 4;
% dotSep = [4,5,6];
% dotSep = linspace(4,6,5);
% Number of canonical orbitals to include
nItinOrbs = 6;
% Dot potential eccentricity (defined as: ly/lx)
eccentricity = 1;
% eccentricity = linspace(0.5,1.5,21);
% Dot bias (bias applied to right dot)
% bias = 0;
bias = linspace(0,0.3,51);
% bias = [0,0.1];
% Tunnel gate potential height (wrt to Vmin)
% Vtun = [1,3,5];
Vtun = 2;
% Vtun = 0;
% Tunnel gate width
lxtun = 2.5/2;

nSweeps = length(nItinOrbs)*length(lxi)*length(dotSep)*length(eccentricity)*...
    length(lxtun)*length(bias)*length(Vtun);

preLoadedCMEs = 1;
simCount = 0;
startTime = clock;
for aa = 1:length(dotSep)
for bb = 1:length(lxi)
for cc = 1:length(nItinOrbs)
for dd = 1:length(eccentricity)
for ee = 1:length(bias)
for ff = 1:length(Vtun)
for gg = 1:length(lxtun)
        simCount = simCount + 1;
        fprintf(1,'*************************\n');
        fprintf(1,'* Simulation: %04d/%04d *\n',simCount,nSweeps);
        fprintf(1,'*************************\n');
        fprintf(1,'* Dot separation = %d    *\n', dotSep(aa));
        fprintf(1,'* Dot width = %.3f     *\n', lxi(bb));
        fprintf(1,'* N canonical orbs = %d *\n', nItinOrbs(cc));
        fprintf(1,'* Dot ecc = %.3f       *\n', eccentricity(dd));
        fprintf(1,'* Dot bias = %.3f     *\n', bias(ee));
        fprintf(1,'* Vtun height = %.3f   *\n', Vtun(ff));
        fprintf(1,'* l_xtun/l_x = %.3f    *\n',lxtun(gg));
        fprintf(1,'*************************\n');
        
        if simCount ~= 1
            % Estimate the remaining runtime
            elapsedTime = etime(clock,startTime);

            avgTimePerLoop = elapsedTime/(simCount-1);
            avgTimeRemaining = avgTimePerLoop*(nSweeps - simCount+1);
            if avgTimeRemaining > 3600
                fprintf(1,'\nTime remaining: %4.1f hours\n\n',avgTimeRemaining/3600);
            elseif avgTimeRemaining > 60
                fprintf(1,'\nTime remaining: %4.1f minutes\n\n',avgTimeRemaining/60);
            else
                fprintf(1,'\nTime remaining: %4.1f seconds\n\n',avgTimeRemaining);
            end
        end
        
        % Start the parpool if it's already not started
%         poolobj = gcp('nocreate'); % If no pool, do not create new one.
%         if isempty(poolobj)
%             parpool('local');
%         end
        

        sparams.dotLocations = [-dotSep(aa)/2,0;dotSep(aa)/2,0];
        [sparams.nDots,~] = size(sparams.dotLocations);
        
        omega = sqrt(Vimin/lxi(bb)^2);

        % Fill in grid parameters
        gparams.ngridx = 200;
        gparams.ngridy = 200;
        gparams.xx = linspace(-8,8,gparams.ngridx)*effaB;
        gparams.yy = linspace(-8,8,gparams.ngridy)*effaB;
        [gparams.XX,gparams.YY] = meshgrid(gparams.xx,gparams.yy);
        currVV = zeros(gparams.ngridy,gparams.ngridx);

        wx = lxi(bb);
        wy = lxi(bb)*eccentricity(dd);
        dot1Pot = (-Vimin-bias(ee)/2)*exp(-((gparams.XX - sparams.dotLocations(1,1)).^2/wx^2) -...
                (gparams.YY - sparams.dotLocations(1,2)).^2/wy^2);
            
        dot2Pot = (-Vimin+bias(ee)/2)*exp(-((gparams.XX - sparams.dotLocations(2,1)).^2/wx^2) -...
                (gparams.YY - sparams.dotLocations(2,2)).^2/wy^2);
            
%         wxtun = wx*lxtunFrac(gg);
        wxtun = lxtun(gg);
        wytun = wy;
        tunPot = Vtun(ff)*exp(-(gparams.XX).^2/wxtun^2 +...
                -(gparams.YY).^2/wytun^2);
            
        gparams.VV = dot1Pot + dot2Pot + tunPot;
        
%         plotMeshgrid(gparams, gparams.VV);
%         continue;
        
        sparams.maxOriginHOsX = 14;
        sparams.maxOriginHOsY = 14;
        sparams.nItinerantOrbitals = nItinOrbs(cc);
        sparams.numElectrons = 2;
    %     sparams.spinSubspaces = 'all';
        sparams.spinSubspaces = [2];
        sparams.nOutputtedEnergies = 4;

        sparams.nOriginHOs = sparams.maxOriginHOsX*sparams.maxOriginHOsY;
        
        tic;
        optOmegaFlag = 1;
        debugFlag = 0;
        if simCount == 1
            if preLoadedCMEs
                [eVecs, ens] = calculateManyBodySpectra_2Bases(sparams, gparams, optOmegaFlag, debugFlag, CMEs_lib);
            else
                [eVecs, ens, CMEs_lib] = calculateManyBodySpectra_2Bases(sparams, gparams, optOmegaFlag, debugFlag);
            end
            exchange = zeros(length(dotSep),length(lxi),...
                length(nItinOrbs),length(eccentricity),length(lxtun),...
                length(bias),length(Vtun));
            energies = zeros(sparams.nOutputtedEnergies,length(dotSep),length(lxi),...
                length(nItinOrbs),length(eccentricity),length(lxtun),...
                length(bias),length(Vtun));
            eVectors = cell(length(dotSep),length(lxi),length(nItinOrbs),...
                length(eccentricity),length(lxtun),length(bias),length(Vtun));
        else
            [eVecs, ens,~,SEens] = calculateManyBodySpectra_2Bases(sparams, gparams, optOmegaFlag, debugFlag, CMEs_lib);
        end
        toc;
        
        exchange(aa,bb,cc,dd,gg,ee,ff) = ens(2,2) - ens(1,1);
        energies(:,aa,bb,cc,dd,gg,ee,ff) = diag(ens);
        eVectors{aa,bb,cc,dd,gg,ee,ff} = eVecs;
end
end
end
end
end
end
end
%%
save dat_biasVSlxi_Vtun=2_Vmin=10_nItin=6.mat sparams exchange energies eVectors lxi Vimin dotSep nItinOrbs eccentricity bias lxtunFrac Vtun

%%
normJFlag = 1;
convertUnits = 0;
% analyzeBiasSweeps(sparams, bias, eccentricity, squeeze(exchange), normJFlag);
% analyzeBiasSweeps(sparams, bias, Vtun, squeeze(exchange), normJFlag);
% analyzeBiasSweeps(sparams, bias, dotSep, squeeze(exchange(:,:,:,:,:)), normJFlag);


%%
% Check accuracy of HO convergence







