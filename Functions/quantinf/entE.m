function ent = entE(psi,varargin)

% ENTE  Entropy of entanglement of a bipartite pure state
% requires: TrX.m
% author: Toby Cubitt
% license: GPL2
%
%    entE(psi,dim) returns the entropy of entanglement of bipartite pure
%    state psi, divided into subsystems with dimensions specified by
%    two-component vector dim.
%
%    If only one dimension is specified, a dim(1) x dim(1) system is
%    assumed.  If two dimensions are specified, a dim(1) x dim(2) system
%    is assumed. If three dimensions are specified, a dim(1) x dim (2) x
%    dim(3) system is assumed, and the entanglement is calculated with
%    respect to the partition 13|2. If no DIM argument is supplied, a
%    two-qubit system is assumed.


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

if length(varargin) == 0
  n = 2;
  m = 2;
  l = 1;
else
  dim = varargin{1};
  switch length(dim)
   case 1
    n = dim(1);
    m = n;
    l = 1;
   case 2
    n = dim(1);
    m = dim(2);
    l = 1;
   case 3
    n = dim(1);
    m = dim(2);
    l = dim(3);
  end
end

if length(psi) ~= n*m*l
  error('MATLAB:badopt',['Dimension of psi is not equal to product' ...
		    ' of dimensions of subsystems'])
end

e = real(eig(TrX(psi,2,[n,m,l])));
indx = find((e ~= 0));
log_e = zeros(size(e));
log_e(indx) = log2(e(indx));
ent = -e'*log_e;
