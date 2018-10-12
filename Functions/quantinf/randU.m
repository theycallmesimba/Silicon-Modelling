function U = randU(n)

% RANDU   Generate random unitary matrix
% requires: nothing
% author: Toby Cubitt
% license: GPL2
%
%    U = RANDU(N) generates a random N x N unitary matrix,
%    distributed uniformly according to the Haar measure.


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


X = (randn(n) + i*randn(n))/sqrt(2);
[Q,R] = qr(X);
R = diag(diag(R)./abs(diag(R)));
U = Q*R;
