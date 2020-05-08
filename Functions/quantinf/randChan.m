function E = randChan(dimA,rep,varargin)

% randChan   Random channel
% requires: randV, chanconv
% author: Toby Cubitt
% license: GPL2
%
%    E = randChan(DIM,REP) returns a randomly generated quantum channel
%    with input and output dimensions DIM, in the representation
%    specified by REP.
%
%    E = randChan(DIMA,DIMB,REP) returns a randomly generated quantum
%    channel with input dimension DIMA and output dimension DIMB, in the
%    representation specified by REP.
%
%    E = randChan(DIMA,DIMB,DIME,REP) returns a randomly generated
%    quantum channel with input dimension DIMA, output dimension DIMB and
%    environment dimension (= "Kraus rank") DIME, in the representation
%    specified by REP.
%
%    REP must be one of:
%      'kraus' - a cell array containing a set of Kraus operators
%      'choi'  - a Choi-Jamiolkowski state representative density matrix
%      'choi2' - as for 'choi', but don't apply local filtering operation
%      'linop' - a matrix representation of the linear operator
%      'isom'  - an isometry matrix
%      'state' - a purification of the Choi-Jamiolkowski state
%
%    See also applychan, chanconv and chan2R.


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


if isempty(varargin)
  dimB = dimA;
  dimE = dimA*dimB;
else
  dimB = rep;
  if length(varargin) == 1
    rep = varargin{1};
    dimE = dimA*dimB;
  else
    dimE = varargin{1};
    rep = varargin{2};
  end
end

% generate random isometry, then convert it to required representation
E = chanconv(randV(dimB*dimE,dimA),'isom',rep,[dimA,dimB]);
