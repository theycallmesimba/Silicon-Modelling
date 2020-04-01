function R = pauli(M)
% PAULI  convert matrix to Hilbert-Schmidt (or Pauli) basis
% requires: nothing
% author: Toby Cubitt
% license: GPL2
%
%    R = PAULI(M) converts a 4x4 matrix from the standard basis to
%    the Hilbert-Schmidt basis. So
%
%      M = sum_i,j R(i,j)*kron(s_i,s_j)/4
%
%    where s_i are pauli matrices (s_0 is the 2x2 identity).


%% Copyright (C) 2004-2010 Toby Cubitt
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


sx = [0 1; 1 0];
sy = [0 -i; i 0];
sz = [1 0; 0 -1];
id = eye(2);

s(:,:,1) = id;
s(:,:,2) = sx;
s(:,:,3) = sy;
s(:,:,4) = sz;

for j = 1:4
  for k = 1:4
    R(j,k) = trace(kron(s(:,:,j),s(:,:,k))*M);
  end
end

