function [ full2DLap ] = make2DSELap(sparams, X, Y, V )
%MAKE2DSELAP Forms the 2DLaplacian for the Schrodinger equation using a 5
%point stencil approach
%   X, Y, V are given in a meshgrid form
%   Use http://www4.ncsu.edu/~zhilin/TEACHING/MA402/chapt5.pdf as reference
%   We assume SI units for the laplacian and a natural ordering for the 2D
%   grid.

    % First make the potential term
    [ny,nx] = size(X);
    PE2DLap = speye(nx*ny,nx*ny);
    
    for ii = 1:nx
        for jj = 1:ny
            % Ignore the warning about it being slow as we initialized the
            % matrix as a sparse identity.  Since we are only changing
            % diagonal elements we don't need to worry about extra overhead
            PE2DLap(ii + nx*(jj-1),ii + nx*(jj-1)) = V(jj,ii);
        end
    end
    
    % Now the KE operator
    dx = X(1,2) - X(1,1);
    dy = Y(2,1) - Y(1,1);
    B = sparse(-((2/(dx*dx)) + (2/(dy*dy)))*eye(nx));
    for ii = 1:nx
        if ii > 1
            B(ii-1,ii) = 1/(dx*dx);
            B(ii,ii-1) = 1/(dx*dx);
        end
    end

    I = sparse(1/(dy*dy)*eye(nx));
    Z = sparse(nx,nx);

    KE2DLap = [];
    for ii = 1:ny
        temp = [];
        for jj = 1:ny
            if ii == jj
                temp = [temp B];
            elseif (jj-1 == ii) || (jj+1 == ii)
                temp = [temp I];
            else
                temp = [temp Z];
            end
        end
        KE2DLap = [KE2DLap;temp];
    end
    
    % Multiply by coefficients for SE equation
    if strcmp(sparams.unitsType,'Rydberg')
        KE2DLap = -KE2DLap;
    else
        KE2DLap = -sparams.hbar^2/(2*sparams.me)*KE2DLap;
    end

    full2DLap = PE2DLap + KE2DLap;
end

