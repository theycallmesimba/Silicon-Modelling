function sparams = solveFor2DLocalizedHOs( sparams, X, Y )
%SOLVEFOR2DLOCALIZEDHOS Summary of this function goes here
%   Detailed explanation goes here

    % The eigenfunctions for the fitted LOHO potentials are simply the
    % product of the LOHO wavefunctions along the x and y axes.  So we can use
    % a simple 1D Schrodinger Equation solver to find the wavefunctions for
    % this 2D potential.

    % Solve for the normalized localized 1 dimensional Harmonic orbitals
    sparams.localXHOs(sparams.nDots*sparams.nLocalOrbitals) = oneDimLOHO;
    sparams.localYHOs(sparams.nDots*sparams.nLocalOrbitals) = oneDimLOHO;
    sparams.localHOs(sparams.nSingleOrbitals) = twoDimLOHO;

    for ii = 1:sparams.nDots
        omega = sparams.fittedPotentialParameters(ii,1);
        ax = sparams.fittedPotentialParameters(ii,3);
        ay = sparams.fittedPotentialParameters(ii,4);
        alpha = sqrt(sparams.me*omega/sparams.hbar);
        for jj = 0:(sparams.nLocalOrbitals-1)
            wfXTemp = 1/sqrt(2^jj*factorial(jj))*(alpha^2/pi)^(1/4)*...
                exp(-alpha^2/2*(X(1,:) - ax).^2).*hermiteH(jj,alpha*(X(1,:) - ax));
            wfYTemp = 1/sqrt(2^jj*factorial(jj))*(alpha^2/pi)^(1/4)*...
                exp(-alpha^2/2*(Y(:,1) - ay).^2).*hermiteH(jj,alpha*(Y(:,1) - ay));
            
            eXTemp = sparams.hbar*omega*(jj + 1/2);
            eYTemp = sparams.hbar*omega*(jj + 1/2);
            
            sparams.localXHOs(jj+1 + sparams.nLocalOrbitals*(ii-1)).initialize(wfXTemp,jj,eXTemp,ii);
            sparams.localYHOs(jj+1 + sparams.nLocalOrbitals*(ii-1)).initialize(wfYTemp,jj,eYTemp,ii);
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

                currWFX = sparams.localXHOs(kk+1 + sparams.nLocalOrbitals*(ii-1)).wavefunction;
                currWFY = sparams.localYHOs(jj+1 + sparams.nLocalOrbitals*(ii-1)).wavefunction;
                currEX = sparams.localXHOs(kk+1 + sparams.nLocalOrbitals*(ii-1)).energy;
                currEY = sparams.localYHOs(jj+1 + sparams.nLocalOrbitals*(ii-1)).energy;
                currWF = currWFY*currWFX;

                sparams.localHOs(ll).initialize(currWF,convertMGtoNO(currWF),kk,jj,currEX + currEY,ii);
            end
        end
    end
end
