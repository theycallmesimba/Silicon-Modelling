function r = flip(p)

% FLIP  Spin flip operation
% requires: nothing
% author: Toby Cubitt
% license: GPL2
%
%    FLIP(RHO) carries out the spin flip operation on a 2-qubit
%    density matrix, used in the definition of concurrence:
%
%      flip(rho) = kron(sy,sy)*conj(rho)*kron(sy,sy)


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


S = [0,0,0,-1;0,0,1,0;0,1,0,0;-1,0,0,0];

switch any(size(p) == 1)
 case true
  r = S*conj(p);
 case false
  r = S*conj(p)*S;
end

