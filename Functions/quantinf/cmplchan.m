function S = cmplchan(E,varargin)

% CMPLCHAN   Construct complementary channel
% requires: chanconv, sysexchange
% author: Toby Cubitt
% license: GPL2
%

%    S = cmplchan(E) constructs the complementary channel of quantum
%    channel E, determining the input representation automatically. The
%    complementary channel is output in the same representation as the
%    input.
%
%    S = cmplchan(E,DIM) specifies the input and output dimensions DIMA
%    and DIMB of the channel in the vector DIM = [DIMA, DIMB].
%
%    S = cmplchan(E,FROM) constructs the complementary channel of quantum
%    channel E, given in representation FROM.
%
%    S = cmplchan(E,FROM,TO) constructs the TO representation of the
%    complementary channel of quantum channel E, given in representation
%    FROM.
%
%    S = cmplchan(E,FROM,DIM) and S = cmplchan(E,FROM,TO,DIM) explicitly
%    specify the input and output dimensions DIMA and DIMB of the channel
%    in the vector DIM = [DIMA, DIMB].
%
%    FROM and TO must be one of:
%      'kraus' - a cell array containing a set of Kraus operators
%      'choi'  - a Choi-Jamiolkowski state representative density matrix
%      'linop' - a matrix representation of the linear operator
%      'isom'  - an isometry matrix
%      'state' - a purification of the Choi-Jamiolkowski state
%      'func'  - a function on density matrices which implements a channel
%      'pauli' - a Pauli representation of a qubit channel
%                (FIXME: not implemented)
%
%   See CHANCONV for more details.


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



% sort out arguments
from = ''; to = ''; dim = [];
if length(varargin) == 1
  if ischar(varargin{1});
    from = varargin{1};
  else
    dim = varargin{1};
  end
elseif length(varargin) == 2
  from = varargin{1};
  if ischar(varargin{2});
    to = varargin{2};
  else
    dim = varargin{2};
  end
elseif length(varargin) == 3
  from = varargin{1};
  to = varargin{2};
  dim = varargin{3};
end

if strcmp(to,'');
  to = from;
end


% convert to purification of choi state
if strcmp(from,'');
  if isempty(dim);
    [S,d] = chanconv(E,'state');
  else
    [S,d] = chanconv(E,'state',dim);
  end
else
  if isempty(dim);
    [S,d] = chanconv(E,from,'state');
  else
    [S,d] = chanconv(E,from,'state',dim);
  end
end

% switch output and environment, and convert back
S = sysexchange(S,[2,3],d);
S = chanconv(S,'state',to,d([1,3]));
