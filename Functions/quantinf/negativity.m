function neg = negativity(p,dim)

% NEGATIVITY   Negativity of density matrix
% requires: nothing
% author: Toby Cubitt
% license: GPL2
%
%    NEGATIVITY(RHO,DIM) returns the negativity of density
%    matrix RHO which has subsystems with dimentions specified
%    by vector DIM.
%
%    Negativity is a measure of entanglemend defined as the
%    sum of the negative eigenvalues of the partial
%    transpose of RHO.
%
%    If no subsystem dimensions are supplied, i.e. DIM=[], a 2 x 2
%    bipartite system is assumed.
%
%    If only one is specified,i.e. DIM=[dim1], a dim(1) x dim(1)
%    bipartite system is assumed.
%
%    If three dimensions are specified, i.e. DIM=[dim1,dim2,dim3],
%    a dim1 x dim2 x dim3 system is assumed and the negativity is
%    calculated for the bipartite splitting
%    sys1 + sys3 | sys2.


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

e = real(eig(Tx(p,2,[n,m,l])));
neg = -2*e'*(e<0);

if neg < 0
  neg = 0;
end
