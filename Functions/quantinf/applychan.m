function rho = applychan(E,psi,varargin)

% applychan   Apply quantum channel E to state PSI
% requires: TrX, chanconv, normalise
% author: Toby Cubitt
% license: GPL2
%
%    RHO = applychan(E,PSI) returns the state (density matrix) RHO
%    obtained by applying the quantum channel E to state PSI (a state
%    vector or density matrix). The representation of the channel is
%    inferred automatically.
%
%    RHO = applychan(E,PSI,REP) returns the state (density matrix) RHO
%    obtained by applying the quantum channel E to state PSI (a state
%    vector or density matrix). The representation of the channel E is
%    determined by REP.
%
%    RHO = applychan(E,PSI,DIM) returns the state (density matrix) RHO
%    obtained by applying the quantum channel E to state PSI (a state
%    vector or density matrix). The input and output dimensions DIMA and
%    DIMB of the channel are specified by the vector DIM = [DIMA, DIMB].
%
%    RHO = applychan(E,PSI,REP,DIM) returns the state (density matrix)
%    RHO obtained by applying the quantum channel E to state PSI (a state
%    vector or density matrix). The representation of the channel E is
%    determined by REP, and the input and output dimensions DIMA and DIMB
%    of the channel are specified by the vector DIM = [DIMA, DIMB].
%
%
%    REP must be one of:
%      'kraus' - a cell array containing a set of Kraus operators
%      'choi'  - a Choi-Jamiolkowski state representative density matrix
%      'choi2' - as for 'choi', but don't apply local filtering operation
%      'linop' - a matrix representation of the linear operator
%      'isom'  - an isometry matrix
%      'state' - a purification of the Choi-Jamiolkowski state
%      'func'  - a function on density matrices which implements a channel
%
%    The entries of a Choi-Jamiolkowski state representative should be
%    ordered so that the input system comes first, the output
%    second. I.e. it should be the state resulting from applying the
%    channel to the *second* half of an entangled state. Non-standard
%    Choi-Jamiolkowski states are allowed as input representations;
%    there is no requirement that the reduced state on the first
%    subsystem be the identity (i.e. it can result from applying the
%    channel to the first half of any full-Schmidt-rank entangled
%    state, which need not necessarily be maximally entangled).
%
%    If the 'choi2' option is specified, the local filtering and
%    normalisation operations required to transform to the correct
%    output state when using a non-standard Choi-Jamiolkowski state
%    are *not* applied. This is only relevant if a non-standard
%    Choi-Jamiolkowski state representative is supplied, i.e. its
%    reduced state on the second subsystem is *not* identity. In this
%    case, RHO is a "twisted" form of the channel output.
%
%    The entries of a state vector representation of a channel should be
%    ordered so that the subsystems occur in the following order: input
%    system, output system, environment.
%
%    It is possible for the same matrix to simultaneously be a valid
%    Choi-Jamiolkowski state representative, a valid matrix
%    representation of the linear operator, and even a valid isometry.
%    Therefore the input representation can not always be inferred
%    unambiguously, and a reasonable guess must be made. The input
%    representation is assumed to be a Choi-Jamiolkowski state
%    representative if E is an (unnormalized) density matrix (i.e. is
%    Hermitian and positive, to within a tolerance of 10E-10). Otherwise
%    it is assumed to be a matrix representation of the linear
%    operator. However, to be absolutely sure the correct representation
%    is used, you should specify it explicitly.
%
%
%    See also randChan, chanconv and chan2R.


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


tol = 1E-10;   % FIXME: should be an optional argument


%----- initialization -----

% sort out function arguments
if length(varargin) == 2
  rep = varargin{1};
  dim = varargin{2};
elseif length(varargin) == 1
  if  ischar(varargin{1});
    rep = varargin{1};
  else
    dim = varargin{1};
  end
end


% if input representation is not explicitly specified, infer it from channel
if ~exist('rep') || strcmp(rep,'')
  if iscell(E)
    rep = 'kraus';
  elseif isa(E, 'function_handle')
    rep = 'func';
  elseif isnumeric(E)
    if size(E,2) == 1
      rep = 'state';
    elseif size(E,1) ~= size(E,2)
      rep = 'isom';
    elseif norm(E - E') < 10E-10 && all(real(eig(E)) >= 0)
      warning(['Can not unambiguously infer input representation; ',...
               'guessing Choi-Jamiolkowski state representative']);
      rep = 'choi';
    else
      warning(['Can not unambiguously infer input representation; ',...
               'guessing linear operator matrix representation']);
      rep = 'linop';
    end
  else
    error('E is not a valid representation of a channel');
  end
end

% check REP argument
if ~(strcmp(rep,'choi') || strcmp(rep,'choi2') || strcmp(rep,'linop') ...
      || strcmp(rep,'kraus') || strcmp(rep,'isom') || strcmp(rep,'state') ...
      || strcmp(rep,'func'))
  error(['REP argument must be one of: ',...
         'choi, choi2, kraus, linop, isom, state or func']);
end


% determine dimensions
if strcmp(rep,'choi') || strcmp(rep,'choi2')
  if exist('dim')
    dimA = dim(1);
    dimB = dim(2);
  else
    dimA = sqrt(size(E,1));  % assume dimA = dimB
    dimB = dimA;
  end
elseif strcmp(rep,'linop')
  if exist('dim')
    dimA = dim(1);
    dimB = dim(2);
  else
    dimA = sqrt(size(E,2));
    dimB = sqrt(size(E,1));
  end
elseif strcmp(rep,'kraus')
  dimA = size(E{1},2);
  dimB = size(E{1},1);
  dimE = length(E);
elseif strcmp(rep,'isom')
  dimA = size(E,2);
  if exist('dim')
    dimB = dim(2);
  else
    dimB = dimA;  % assume dimB = dimA
  end
  dimE = size(E,1)/dimB;
elseif strcmp(rep,'state')
  if exist('dim')
    dimA = dim(1);
    dimB = dim(2);
    dimE = length(E)/(dimA*dimB);
  else
    error('State input representation requires explicit DIM argument')
  end
end


% check consistency of dimensions
if ~strcmp(rep,'func')
  if mod(dimA,1) ~= 0 || mod(dimB,1) ~= 0 || ...
        (exist('dimE') && mod(dimE,1) ~= 0)
    error(['Dimensions could not be inferred automatically, ',...
           'or E is not a valid representation of a channel']);
  end
  if exist('dim') && (dimA ~= dim(1) || dimB ~= dim(2))
    error(['DIM argument does not match channel dimensions']);
  end
end


% convert an input state vector into a density matrix
if size(psi,1) == 1
  psi = psi';
end
if any(size(psi) == 1)
  psi = psi*psi';
end




%----- algorithms -----

if strcmp(rep,'kraus')
  rho = zeros(dimB);
  for k = 1:dimE
    rho = rho + E{k}*psi*E{k}';
  end

elseif strcmp(rep,'func')
  rho = E(psi);

elseif strcmp(rep,'choi')
  % convert from choi to choi before applying, to ensure it is a
  % standard state (i.e. reduced state on second system is identity)
  E = chanconv(E,'choi','choi',[dimA,dimB]);
  rho = TrX(E*kron(psi.',eye(dimB)),1,[dimA,dimB]);

elseif strcmp(rep,'choi2')
%  rho = normalise(TrX(E*kron(psi.',eye(dimB)),1,[dimA,dimB]));
  rho = TrX(E*kron(psi.',eye(dimB)),1,[dimA,dimB]);

elseif strcmp(rep,'linop')
  rho = reshape(E*reshape(psi,dimA^2,1),dimB,dimB);

elseif strcmp(rep,'isom') || strcmp(rep,'state')
  if strcmp(rep,'state')
    % convert state to isometry in order to apply it
    E = chanconv(E,'state','isom',[dimA,dimB]);
  end
  rho = TrX(E*psi*E',2,[dimB,dimE]);

end
