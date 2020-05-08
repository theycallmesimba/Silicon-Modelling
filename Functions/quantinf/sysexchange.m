function q = sysexchange(p,sys,dim)

% SYSEXCHANGE  exchange order of two subsystems in a multipartite state
% requires: syspermute
% author: Toby Cubitt
% license: GPL2
%
%    Q = SYSEXCHANGE(P,SYS,DIM) exchanges the order of the two
%    subsystems of state P (a state vector or density matrix)
%    specified by SYS (a vector of two integers). DIM is a vector
%    specifying the dimensions of the component subsystems.


%% Copyright (C) 2010 Toby Cubitt
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

perm = [1:length(dim)];
perm(sys(1)) = sys(2);
perm(sys(2)) = sys(1);

q = syspermute(p,perm,dim);
