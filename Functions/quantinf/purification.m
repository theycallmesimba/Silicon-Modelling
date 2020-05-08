function varargout = purification(rho,varargin)

% PURIFICATION
% requires: nothing
% author: Toby Cubitt
% license: GPL2
%
%    PSI = PURIFICATION(RHO) returns a minimal purification of density matrix
%          RHO. Eigenvalues of RHO smaller than a tolerance of
%
%              tol = length(RHO) * sigma(1) * eps
%
%          are taken to be zero, where sigma(1) is the largest singular
%          value of RHO and eps is the machine precision.
%
%    PSI = PURIFICATION(RHO,TOL) sets the tolerance to TOL.
%
%    [PSI DIM] = PURIFICATION(RHO) returns the dimension of the auxilliary
%    purifying system as DIM.


%% Copyright (C) 2004-2009 Toby Cubitt
%%
%% This program is free software; you can redistribute it and/or
%% modify it under the terms of the GNU General Public License
%% as published by the Free Software Foundation; either version 2
%% of the License, or (at your option) any later version.
%%
%% This program is distributed in the hope that it will be useful,
%% but WITHOUT ANY WARRANTY; without even the implied warranty of
%% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%% GNU General Public License for more details.
%%
%% You should have received a copy of the GNU General Public License
%% along with this program; if not, write to the Free Software
%% Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
%% MA 02110-1301, USA.


% sort out input arguments
dim = length(rho);
if isempty(varargin)
  sigma = svd(rho);
  tol = dim*sigma(1)*eps;
else
  tol = varargin(1);
end


% calculate purification by forming direct sum of weighted eigenvectors
[U,D] = eig(rho);
[D,idx] = sort(diag(D),'descend');
U = U(:,idx);
dimE = sum(D > tol);
psi = zeros(dim*dimE,1);
for j = 1:dimE
  psi([j:dimE:dim*dimE]) = sqrt(D(j))*U(:,j);
end


% sort out output arguments
varargout{1} = psi;
if nargout > 1
  varargout{2} = dimE;
end
