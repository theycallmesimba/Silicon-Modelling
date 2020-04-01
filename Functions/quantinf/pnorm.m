
function n = pnorm(M,p)

% PNORM   Schatten p-norm
% requires: nothing
% author: Toby Cubitt
% license: GPL2
%
%    N = PNORM(M,p) computes the Schatten p-norm of M.
%    p = 1 is the trace norm.
%    p = 2 is the Frobenius norm (note that norm(M,2) is more efficient).
%    p = 'inf' is the spectral norm.
%    p < 1 are pseudo-norms.
%    p = 0 is the rank.


%% Copyright (C) 2012 Toby Cubitt
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

if p == 'inf'
   n = max(svd(M));
elseif p == 0
  n = rank(M);
else
  n = (sum(svd(M).^p))^(1/p);
end
