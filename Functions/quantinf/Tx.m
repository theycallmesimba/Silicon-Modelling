function x = Tx(rho,sys,dim)

% Tx   Partial transpose
% requires: nothing
% author: Toby Cubitt
% license: GPL2
%
%    X = Tx(RHO,SYS,DIM) gives the partial transpose of a matrix RHO
%    with respect to subsystems SYS, where the subsystem dimensions are
%    specified by vector DIM
%
%    If only one dimension is specified, i.e. DIM=dim1, a dim1 x dim1
%    system is assumed.


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


n = length(dim);
d = size(rho);
perm = [1:2*n];
perm([n+1-sys,2*n+1-sys]) = perm([2*n+1-sys,n+1-sys]);
x = reshape(permute(reshape(rho,[dim(end:-1:1),dim(end:-1:1)]),perm),d);
