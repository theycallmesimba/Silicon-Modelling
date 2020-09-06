function LOHOs = solveFor2DLocalizedHOs( sparams, gridparams )
%SOLVEFOR2DLOCALIZEDHOS Summary of this function goes here
%   Detailed explanation goes here

    XX = gridparams.XX;
    YY = gridparams.YY; 
    
    % The eigenfunctions for the fitted LOHO potentials are simply the
    % product of the LOHO wavefunctions along the x and y axes.

    % Solve for the normalized localized 1 dimensional Harmonic orbitals
    XLOHOs(sparams.nDots*sparams.nLocalOrbitals) = oneDimLOHO;
    YLOHOs(sparams.nDots*sparams.nLocalOrbitals) = oneDimLOHO;
    LOHOs(sparams.nSingleOrbitals) = twoDimLOHO;

    for ii = 1:sparams.nDots
        sparams.omegaLOHO(ii) = sparams.fittedPotentialParameters(ii,1);
        ax = sparams.fittedPotentialParameters(ii,3);
        ay = sparams.fittedPotentialParameters(ii,4);
        alpha = sqrt(sparams.me*sparams.omegaLOHO(ii)/sparams.hbar);
        for jj = 0:(sparams.nLocalOrbitals-1)
            wfXTemp = 1/sqrt(2^jj*factorial(jj))*(alpha^2/pi)^(1/4)*...
                exp(-alpha^2/2*(XX(1,:) - ax).^2).*hermiteH(jj,alpha*(XX(1,:) - ax));
            wfYTemp = 1/sqrt(2^jj*factorial(jj))*(alpha^2/pi)^(1/4)*...
                exp(-alpha^2/2*(YY(:,1) - ay).^2).*hermiteH(jj,alpha*(YY(:,1) - ay));
            
            eXTemp = sparams.hbar*sparams.omegaLOHO(ii)*(jj + 1/2);
            eYTemp = sparams.hbar*sparams.omegaLOHO(ii)*(jj + 1/2);
            
            XLOHOs(jj+1 + sparams.nLocalOrbitals*(ii-1)).initialize(wfXTemp,jj,eXTemp,ii);
            YLOHOs(jj+1 + sparams.nLocalOrbitals*(ii-1)).initialize(wfYTemp,jj,eYTemp,ii);
        end
    end

    % Construct the local phi(x,y) HO wavefunctions
    ll = 0;
    for ii = 1:sparams.nDots
        for jj = 0:(sparams.nLocalOrbitals-1)
            for kk = 0:(sparams.nLocalOrbitals-1)
                
                if jj + kk >= sparams.nLocalOrbitals
                    continue;
                end
                ll = ll + 1;

                currWFX = XLOHOs(kk+1 + sparams.nLocalOrbitals*(ii-1)).wavefunction;
                currWFY = YLOHOs(jj+1 + sparams.nLocalOrbitals*(ii-1)).wavefunction;
                currEX = XLOHOs(kk+1 + sparams.nLocalOrbitals*(ii-1)).energy;
                currEY = YLOHOs(jj+1 + sparams.nLocalOrbitals*(ii-1)).energy;
                currWF = currWFY*currWFX;

                LOHOs(ll).initialize(currWF,convertMGtoNO(currWF),kk,jj,currEX + currEY,ii);
            end
        end
    end
end
