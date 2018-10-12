function q = syspermute(p,perm,dim)

% SYSPERMUTE  permute order of subsystems in a multipartite state
% requires: nothing
% author: Toby Cubitt
% license: GPL2
%
%    Q = SYSPERMUTE(P,PERM,DIM) permutes the order pf the subsystems
%    of state P (a state vector or density matrix) according to
%    permutation PERM. DIM is a vector specifying the dimensions of
%    the component subsystems. PERM should contain the numbers 1
%    through length(DIM) in the order in which the subsystems should
%    appear in the output state Q.


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


% check arguments
n = length(dim);
d = size(p);
if length(perm) ~= n
  error('Number of subsystems in PERM and DIM inconsistent')
end
if sort(perm) ~= [1:n]
  error('PERM is not a valid permutation')
end
if length(p) ~= prod(dim)
  error('Total dimension in DIM does not match state P')
end


if any(size(p)==1)
  % state vector
  perm = n+1-perm([end:-1:1]);
  q = reshape(permute(reshape(p,dim(end:-1:1)),perm),d);

elseif d(1) == d(2)
  % density matrix
  perm = n+1-perm([end:-1:1]);
  perm = [perm,n+perm];
  q = reshape(permute(reshape(p,[dim(end:-1:1),dim(end:-1:1)]),perm),d);
end
