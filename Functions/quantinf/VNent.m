function x = VNent(p)

% VNENT   von-Neumann entropy
% requires: nothing
% author: Toby Cubitt
% license: GPL2
%
%     VNENT(RHO) gives the von-Neumann entropy of density matrix RHO
%     (note: this is NOT the entropy of the reduced density matrix)


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


e = eig(p);
x = -e'*log2(e+(e==0));

