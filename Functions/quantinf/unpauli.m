function M = unpauli(R)

% UNPAULI  convert matrix from Hilbert-Schmidt (or Pauli) basis
% requires: nothing
% author: Toby Cubitt
% license: GPL2
%
%    M = UNPAULI(R) converts a 4x4 matrix from the Hilbert-Schmidt
%    basis back to the standard basis. So
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

M = zeros(4);
for j = 1:4
  for k = 1:4
    M = M+R(j,k)*kron(s(:,:,j),s(:,:,k))/4;
  end
end

