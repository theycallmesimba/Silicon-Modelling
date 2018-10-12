function R = chan2R(S)

% chan2R   Convert the Jamiolkowski representation of a qubit channel to
%          the corresponding R representation
% requires: pauli
% author: Toby Cubitt
% license: GPL2
%
%    R = chan2R(S) returns the R representation of a qubit channel S, where S
%    is the channel's Jamiolkowski representation.
%
%    See also chanconv, randChan and applychan.


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


R = pauli(Tx(S,2,[2,2]));
