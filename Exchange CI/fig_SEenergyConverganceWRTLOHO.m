nHOOrbsStart = 5;
nHOOrbsEnd = 8;
nHOOrbs = nHOOrbsStart:nHOOrbsEnd;

lxi = 2.3;
Vimin = 10;
dotSep = 5;
eccentricity = 1.2;
bias = 0.4;
Vtun = 2;
nItinOrbs = 10;
lxtunFrac = 1/2;

manyBodyEns = zeros(length(nHOOrbs),nItinOrbs);
manyBodyEnsOptOmega = zeros(length(nHOOrbs),nItinOrbs);
count = 1;
for jj = [0,1]
    for qq = 1:length(nHOOrbs)        
        sparams.dotLocations = [-dotSep/2,0;dotSep/2,0];
        [sparams.nDots,~] = size(sparams.dotLocations);
        
        omega = sqrt(Vimin/lxi^2);
        if count == 1
            % Fill in grid parameters
            gparams.ngridx = 200;
            gparams.ngridy = 200;
            gparams.xx = linspace(-8,8,gparams.ngridx)*effaB;
            gparams.yy = linspace(-8,8,gparams.ngridy)*effaB;
            [gparams.XX,gparams.YY] = meshgrid(gparams.xx,gparams.yy);
            currVV = zeros(gparams.ngridy,gparams.ngridx);

            wx = lxi;
            wy = lxi*eccentricity;
            dot1Pot = -Vimin*exp(-((gparams.XX - sparams.dotLocations(1,1)).^2/wx^2) -...
                    (gparams.YY - sparams.dotLocations(1,2)).^2/wy^2);

            dot2Pot = (-Vimin+bias)*exp(-((gparams.XX - sparams.dotLocations(2,1)).^2/wx^2) -...
                    (gparams.YY - sparams.dotLocations(2,2)).^2/wy^2);

            wxtun = wx*lxtunFrac;
            wytun = wy;
            tunPot = Vtun*exp(-(gparams.XX).^2/wxtun^2 +...
                    -(gparams.YY).^2/wytun^2);

            gparams.VV = dot1Pot + dot2Pot + tunPot;
        end
            
        sparams.nItinerantOrbitals = nItinOrbs;
        sparams.maxOriginHOsX = nHOOrbs(qq);
        sparams.maxOriginHOsY = nHOOrbs(qq);
        sparams.nOriginHOs = sparams.maxOriginHOsX*sparams.maxOriginHOsY;
        
        omegaGuess = abs(mean(sparams.fittedPotentialParameters(:,1)));
        if jj == 0
            optOmega = omegaGuess;
        elseif jj == 1
            fprintf(1,'Optimizing origin harmonic orbital omega...  \n');
            optOmega = optimizeOmega(sparams,gparams,omegaGuess);
            fprintf(1,'Done!\n\n');
        end
        
        fprintf(1,'Finding 2D harmonic orbitals at origin...  \n');
        % Create a new basis of HOs centered at the origin of our dots
        [originHOs, originOmega] = createOriginHOs(sparams,gparams,optOmega);
        
        nStates = sparams.nOriginHOs;
        basisToUse = originHOs;
        
        if count == 1
            [itinOrbs,itinEns] = findItinerantBasis(sparams, gparams, sparams.nItinerantOrbitals);
        end
        acoeffs = findTMatrixViaInnerProd(gparams, basisToUse, itinOrbs);

        % Now find the itinerant basis energies using acoeffs
        sparams.LCHOEnergies = zeros(1,sparams.nItinerantOrbitals);
        full2DLap = make2DSELap(sparams,gparams);
        for kk = 1:sparams.nItinerantOrbitals
            tempwf = zeros(gparams.ngridy*gparams.ngridx,1);

            for ll = 1:nStates
                tempwf = tempwf + acoeffs(kk,ll)*originHOs(ll).wavefunctionNO;
            end

            sparams.LCHOEnergies(kk) = getInnerProduct2D(itinOrbs(kk).wavefunctionMG,...
                convertNOtoMG(full2DLap*tempwf,gparams.ngridx,gparams.ngridy), gparams.XX, gparams.YY);
        end

        % Truncate A in accordance with how many itinerant orbitals we want
        acoeffs = acoeffs(1:sparams.nItinerantOrbitals,:);
        sparams.LCHOEnergies = sparams.LCHOEnergies(1:sparams.nItinerantOrbitals);
        
        ens = sort(sparams.LCHOEnergies);
        if jj == 0
            manyBodyEns(qq,:) = ens(1:sparams.nItinerantOrbitals);
        elseif jj == 1
            manyBodyEnsOptOmega(qq,:) = ens(1:sparams.nItinerantOrbitals);
        end
        count = 2;
        
%         clearvars sparams;
    end
end

fprintf(1,'Done with convergence simulation.\n');
%%
cmColors = viridis(nItinOrbs);
figure('Color','white');
hold on;
fullItinEns = repmat(itinEns,1,length(nHOOrbs))';
perDiffNaiveOmega = abs((manyBodyEnsOptOmega - fullItinEns)./fullItinEns);
perDiffOptOmega = abs((manyBodyEns - fullItinEns)./fullItinEns);
% hG3Lines = plot([nHOOrbsStart,nHOOrbsEnd],repmat(itinEns,1,2),'Color','k','Linestyle','--','Linewidth',2);
% hG1Lines = plot(nHOOrbsStart:nHOOrbsEnd,manyBodyEns(1:length(nHOOrbs),:),'Color',cm(20,:),'linewidth',2);
% hG2Lines = plot(nHOOrbsStart:nHOOrbsEnd,manyBodyEnsOptOmega(1:length(nHOOrbs),:),'Color',cm(110,:),'linewidth',2);
for ii = 1:nItinOrbs
    plot(nHOOrbsStart:nHOOrbsEnd,perDiffNaiveOmega(:,ii),'Color',cmColors(ii,:),'Linestyle','--','Linewidth',2);
    plot(nHOOrbsStart:nHOOrbsEnd,perDiffOptOmega(:,ii),'Color',cmColors(ii,:),'Linestyle','-','Linewidth',2);
end
% hG1Lines = plot(nHOOrbsStart:nHOOrbsEnd,perDiffNaiveOmega,'Color',cm(20,:),'linewidth',2);
% hG2Lines = plot(nHOOrbsStart:nHOOrbsEnd,perDiffOptOmega,'Color',cm(110,:),'linewidth',2);
set(gca,'TickLabelInterpreter','latex','Fontsize',14,'YScale','log');
xlabel('$M_x, M_y = \sqrt{M}$','Interpreter','latex','Fontsize',20);
ylabel('$\epsilon''_j$ [Ry*]','Interpreter','latex','Fontsize',20);
ylabel('$|\frac{\epsilon_i'' - \epsilon_i}{\epsilon_i}|$','Interpreter','latex','Fontsize',20);
xlim([nHOOrbsStart,nHOOrbsEnd]);
% xlim([9,nOrbsEnd]);
% ylim([-24.5,-23.9]);
% title('Energy convergence: 10 Canons','Interpreter','latex','Fontsize',20);
% title('Energy convergence','Interpreter','latex','Fontsize',20);
% hG1Group = hggroup;
% hG2Group = hggroup;
% hG3Group = hggroup;
% set(hG3Lines,'Parent',hG3Group);
% set(hG1Lines,'Parent',hG1Group)
% set(hG2Lines,'Parent',hG2Group)
% Include these hggroups in the legend:
% set(get(get(hG3Group,'Annotation'),'LegendInformation'),...
%     'IconDisplayStyle','on'); 
% set(get(get(hG1Group,'Annotation'),'LegendInformation'),...
%     'IconDisplayStyle','on'); 
% set(get(get(hG2Group,'Annotation'),'LegendInformation'),...
%     'IconDisplayStyle','on'); 
% legend({'without $\omega$ optimization','with $\omega$ optimization','$\epsilon_j$'},'Location','northeast','Interpreter','latex');
legend({'without $\omega$ optimization','with $\omega$ optimization'},'Location','northeast','Interpreter','latex')
legend('boxoff')
% ylim([-8,-1]);
% export_fig 'SEconvWRTloho.png' -m5

%%
nEns = 6;
cm = viridis();
figure('Color','white');
hold on;
fullItinEns = repmat(itinEns,1,length(nHOOrbs))';
perDiffNaiveOmega = abs((manyBodyEnsOptOmega - fullItinEns)./fullItinEns);
perDiffOptOmega = abs((manyBodyEns - fullItinEns)./fullItinEns);
hG3Lines = plot([nHOOrbsStart,nHOOrbsEnd],repmat(itinEns(1:nEns),1,2),'Color','k','Linestyle','--','Linewidth',2);
hG1Lines = plot(nHOOrbsStart:nHOOrbsEnd,manyBodyEns(1:length(nHOOrbs),1:nEns),'Color',cm(20,:),'linewidth',2);
hG2Lines = plot(nHOOrbsStart:nHOOrbsEnd,manyBodyEnsOptOmega(1:length(nHOOrbs),1:nEns),'Color',cm(110,:),'linewidth',2);
set(gca,'TickLabelInterpreter','latex','Fontsize',14,'YScale','log');
xlabel('$M_x, M_y = \sqrt{M}$','Interpreter','latex','Fontsize',20);
ylabel('$\epsilon''_j$ [Ry*]','Interpreter','latex','Fontsize',20);
xlim([nHOOrbsStart,nHOOrbsEnd]);

hG1Group = hggroup;
hG2Group = hggroup;
hG3Group = hggroup;
set(hG3Lines,'Parent',hG3Group);
set(hG1Lines,'Parent',hG1Group)
set(hG2Lines,'Parent',hG2Group)
% Include these hggroups in the legend:
set(get(get(hG3Group,'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','on'); 
set(get(get(hG1Group,'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','on'); 
set(get(get(hG2Group,'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','on'); 
legend({'without $\omega$ optimization','with $\omega$ optimization','$\epsilon_j$'},'Location','northeast','Interpreter','latex');
legend('boxoff')

% export_fig 'SEconvWRTloho.png' -m5





