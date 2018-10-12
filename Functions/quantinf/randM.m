function M = randM(varargin)

% RANDM   Generate random complex matrix
% requires: nothing
% author: Toby Cubitt
% license: GPL2
%
%    M = RANDM(N) generates a random N x N complex matrix.
%
%    M = RANDM(M,N,P,...) generates a random M x N x P x ... complex
%    tensor.


%% Copyright (C) 2004-2011 Toby Cubitt
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


M = 2*(rand(varargin{:}) + i*rand(varargin{:})) - (1+i);
