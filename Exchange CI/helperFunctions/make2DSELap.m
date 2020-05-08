function [ full2DLap ] = make2DSELap(sparams, gparams)
%MAKE2DSELAP Forms the 2DLaplacian for the Schrodinger equation using a 5
%point stencil approach
%   Use http://www4.ncsu.edu/~zhilin/TEACHING/MA402/chapt5.pdf as reference
%   We allow for both SI and Rydberg units for the laplacian and a natural 
%   ordering for the 2D grid.

    % First make the potential term
    nx = gparams.ngridx;
    ny = gparams.ngridy;
%     PE2DLap = speye(nx*ny,nx*ny);
    
    % Create a diagonal matrix with the diagonal elements written out in
    % the natural ordering for the meshgrid.
    PE2DLap = sparse(1:nx*ny,1:nx*ny,convertMGtoNO(gparams.VV));
    
%     for ii = 1:nx
%         for jj = 1:ny
%             % Ignore the warning about it being slow as we initialized the
%             % matrix as a sparse identity.  Since we are only changing
%             % diagonal elements we don't need to worry about extra overhead
%             PE2DLap(ii + nx*(jj-1),ii + nx*(jj-1)) = V(jj,ii);
%         end
%     end
    
    % Now the KE operator
    dx = gparams.XX(1,2) - gparams.XX(1,1);
    dy = gparams.YY(2,1) - gparams.YY(1,1);

    B = speye(nx,nx)*(-(2/(dx*dx) + 2/(dy*dy)));
    % Add +/-1 offdiagonal entries for the 1/dx^2 elements
    B = B + sparse(diag(ones(nx-1,1)/(dx*dx),-1));
    B = B + sparse(diag(ones(nx-1,1)/(dx*dx),1));
    
    % Now create the block diagonal matrix of Bs..
    KE2DLap = kron(speye(ny),B);
    % Now set the off diagonal entries for the 1/dy^2 elements
    KE2DLap = KE2DLap + kron(sparse(diag(ones(ny-1,1),-1)),speye(nx)/(dy*dy));
    KE2DLap = KE2DLap + kron(sparse(diag(ones(ny-1,1),1)),speye(nx)/(dy*dy));
    
%     for ii = 1:nx
%         if ii > 1
%             B(ii-1,ii) = 1/(dx*dx);
%             B(ii,ii-1) = 1/(dx*dx);
%         end
%     end
% 
%     I = sparse(1/(dy*dy)*eye(nx));
%     Z = sparse(nx,nx);
% 
%     KE2DLap = [];
%     for ii = 1:ny
%         temp = [];
%         for jj = 1:ny
%             if ii == jj
%                 temp = [temp B];
%             elseif (jj-1 == ii) || (jj+1 == ii)
%                 temp = [temp I];
%             else
%                 temp = [temp Z];
%             end
%         end
% 
%         KE2DLap = [KE2DLap;temp];
%     end
    
    % Multiply by coefficients for SE equation
    if strcmp(sparams.unitsType,'Rydberg')
        KE2DLap = -KE2DLap;
    else
        KE2DLap = -sparams.hbar^2/(2*sparams.me)*KE2DLap;
    end
    full2DLap = PE2DLap + KE2DLap;        
end