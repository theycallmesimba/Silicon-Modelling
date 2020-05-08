function [ full1DLap ] = make1DSELap( sparams, X, V )
%MAKE2DSELAP Forms the 2DLaplacian for the Schrodinger equation using a 3
%point stencil approach
%   X, V are linear arrays
%`  We assume that X has uniform spacing between points

    % Write potential in matrix form
    nv = length(V);
    PE1DLap = speye(nv,nv);
    for ii = 1:nv
        PE1DLap(ii,ii) = V(ii);
    end
    
%     % Form the matrix for the KE operator (3 point stencil)
%     KE1DLap = speye(nv,nv);
%     for ii = 1:nv
%         KE1DLap(ii,ii) = -2;
%         if ii > 1
%             KE1DLap(ii-1,ii) = 1;
%             KE1DLap(ii,ii-1) = 1;
%         end
%     end
    % Form the matrix for the KE operator (5 point stencil)
    KE1DLap = speye(nv,nv);
    for ii = 1:nv
        KE1DLap(ii,ii) = -30;
        if ii > 1
            KE1DLap(ii-1,ii) = 16;
            KE1DLap(ii,ii-1) = 16;
        end
        if ii > 2
            KE1DLap(ii-2,ii) = -1;
            KE1DLap(ii,ii-2) = -1;
        end
    end
    KE1DLap = KE1DLap/12;
    
    dx = X(2) - X(1);
    KE1DLap = -sparams.hbar^2/(2*sparams.me*dx*dx)*KE1DLap;

    % Now form the whole Laplacian and solve for wave functions
    full1DLap = KE1DLap + PE1DLap;
end

