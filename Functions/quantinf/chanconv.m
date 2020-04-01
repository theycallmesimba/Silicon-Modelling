function varargout = chanconv(E,from,varargin)

% CHANCONV   Convert between representations of a quantum channel
% requires: TrX, schmidt, purification
% author: Toby Cubitt
% license: GPL2
%
%    S = chanconv(E,TO) converts the quantum channel E to the TO
%    representation, determining the input representation automatically.
%
%    S = chanconv(E,FROM,TO) converts the quantum channel E in the FROM
%    representation to the TO representation.
%
%    S = chnconv(E,TO,DIM) converts the quantum channel E to the TO
%    representation, determining the input representation
%    automatically. The input and output dimensions DIMA and DIMB of the
%    channel are specified by the vector DIM = [DIMA, DIMB].
%
%    S = chnconv(E,FROM,TO,DIM) converts the quantum channel E in the
%    FROM representation to the TO representation. The input and output
%    dimensions DIMA and DIMB of the channel are specified by the vector
%    DIM = [DIMA, DIMB].
%
%    [S,dim] = chancon(E, ...) returns the input, output and environment
%    dimensions of the channel in a vector DIM = [DIMA,DIMB,DIME].
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
%    The Kraus operators are the operators A_k acting from the *left* in
%    the Kraus decomposition (i.e. the action of the channel is given by
%    sum_k A_k rho A_k').
%
%    The entries of a Choi-Jamiolkowski state representative matrix
%    should be ordered so that the input system comes first, the output
%    second. I.e. it should be the state resulting from applying the
%    channel to the *second* half of an entangled state. Non-standard
%    Choi-Jamiolkowski states are allowed as input representations; there
%    is no requirement that the reduced state on the first subsystem be
%    the identity (i.e. it can result from applying the channel to the
%    first half of any full-Schmidt-rank entangled state, which need not
%    necessarily be maximally entangled).
%
%    Output Choi-Jamiolkowski state representatives follow the same
%    ordering convention, but the reduced state on the first subsystem is
%    always the identity (i.e. they are always the state resulting from
%    applying the channel to the second half of a canonical maximally
%    entangled state).
%
%    The entries of a state vector representation of a channel should be
%    ordered so that the subsystems occur in the following order: input
%    system, output system, environment. Ouput state vectors follow the
%    same ordering convention.
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
%    operator. However, to be absolutely sure the correct input
%    representation is used, you should specify it explicitly.
%
%    See also chan2R, randChan and applychan.


%% Copyright (C) 2004-2010, 2012 Toby Cubitt
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



%----- initialization -----

% sort out function arguments
if length(varargin) == 2
  to = varargin{1};
  dim = varargin{2};
elseif length(varargin) == 1
  if  ischar(varargin{1});
    to = varargin{1};
  else
    dim = varargin{1};
    to = from;
    from = '';
  end
else
  to = from;
  from = '';
end


% if input representation is not explicitly specified, infer it from channel
if strcmp(from,'')
  if iscell(E)
    from = 'kraus';
  elseif isa(E,'function_handle')
    from = 'func';
  elseif isnumeric(E)
    if size(E,2) == 1
      from = 'state';
    elseif size(E,1) ~= size(E,2)
      from = 'isom';
    elseif norm(E - E') < 10E-10 && all(real(eig(E)) >= 0)
      warning(['Can not unambiguously infer input representation; guessing Choi-Jamiolkowski state representative']);
      from = 'choi';
    else
      warning(['Can not unambiguously infer input representation; guessing linear operator matrix representation']);
      from = 'linop';
    end
  else
    error('E is not a valid representation of a channel');
  end
end

% check FROM and TO representation arguments
if ~(strcmp(to,'choi') || strcmp(to,'linop') ...
     || strcmp(to,'kraus') || strcmp(to,'isom') ...
     || strcmp(to,'state') || strcmp(to,'func')) ...
   || ~(strcmp(from,'choi') || strcmp(from,'linop') ...
     || strcmp(from,'kraus') || strcmp(from,'isom') ...
     || strcmp(from,'state') || strcmp(from,'func'))
  error(['FROM and TO arguments must be one of: ',...
         'choi, kraus, linop, isom, state or func']);
end


% determine dimensions
if strcmp(from,'choi')
  if exist('dim')
    dimA = dim(1);
    dimB = dim(2);
  else
    dimA = sqrt(size(E,1));  % assume dimA = dimB
    dimB = dimA;
  end
elseif strcmp(from,'linop')
  if exist('dim')
    dimA = dim(1);
    dimB = dim(2);
  else
    dimA = sqrt(size(E,2));
    dimB = sqrt(size(E,1));
  end
elseif strcmp(from,'kraus')
  dimA = size(E{1},2);
  dimB = size(E{1},1);
  dimE = length(E);
elseif strcmp(from,'isom')
  dimA = size(E,2);
  if exist('dim')
    dimB = dim(2);
  else
    dimB = dimA;  % assume dimB = dimA
  end
  dimE = size(E,1)/dimB;
elseif strcmp(from,'state')
  if exist('dim');
    dimA = dim(1);
    dimB = dim(2);
    dimE = length(E)/(dimA*dimB);
  else
    error('State input representation requires explicit DIM argument');
  end
elseif strcmp(from,'func')
  if exist('dim');
    dimA = dim(1);
    if length(dim) > 1
      dimB = dim(2);
    end
  else
    error(['Function input representation requires explicit DIM argument']);
  end
end

% check consistency of dimensions
if mod(dimA,1) ~= 0 || (exist('dimB') && mod(dimB,1)) ~= 0 || ...
      (exist('dimE') && mod(dimE,1) ~= 0)
  error(['Dimensions could not be inferred automatically, ', ...
         'or E is not a valid representation of a channel']);
end
if exist('dim') && (dimA ~= dim(1) || (exist('dimB') && dimB ~= dim(2)))
  error(['DIM argument does not match channel dimensions']);
end




%----- algorithms -----

% set ouput to input, in case we're converting from and to the same
% representation
S = E;


%--- kraus ---
if strcmp(from,'kraus')
  if strcmp(to,'func')
    % convert from kraus to func by building a function that calls 'applychan'
    S = @(rho) applychan(S,rho,'kraus',[dimA,dimB]);

  elseif strcmp(to,'isom') || strcmp(to,'state')
    % convert to isometry by tensoring with projectors on computational
    % basis states of the environment; if converting to state, convert to
    % isometry and leave rest till later
    S = zeros(dimB*dimE,dimA);
    for j = 1:dimE
      S([j:dimE:dimB*dimE],:) = E{j};
    end
    from = 'isom';
    E = S;

  elseif strcmp(to,'choi') || strcmp(to,'kraus')
    % convert from kraus to choi by applying the channel to second half
    % of a maximally entangled state; if converting from kraus to kraus,
    % convert to choi and convert back later, to produce a minimal kraus
    % representation
    w = zeros(dimA^2,1);
    w([1:dimA+1:dimA^2]) = 1;
    w = w*w';
    S = zeros(dimA*dimB);
    for j = 1:dimE
      S = S + kron(eye(dimA),E{j})*w*kron(eye(dimA),E{j}');
    end
    from = 'choi';
    E = S;

  elseif strcmp(to,'linop')
    % convert kraus to linop by summing all tensor products of kraus
    % operators with their complex conjugates; if converting to anything
    % else, convert to linop and leave rest till later
    S = zeros(dimB^2,dimA^2);
    for j = 1:length(E)
      S = S + kron(conj(E{j}),E{j});
    end
    from = 'linop';
    E = S;
  end
end



%--- func ---
if strcmp(from,'func') && ~strcmp(to,'func')
  % construct linear operator by applying function to complete operator
  % basis, and leave rest till later
  clear('S');
  for j = 1:dimA
    for k = 1:dimA
      e = zeros(dimA);
      e(k,j) = 1;
      e = E(e);
      S(:,dimA*(j-1)+k) = reshape(e,[prod(size(e)),1]);
    end
  end

  % calculate / check consistency of output dimension
  if ~exist('dimB')
    dimB = sqrt(size(S,1));
    if mod(dimB,1) ~= 0
      error(['function E does not implement a channel']);
    end
  elseif dimB ~= sqrt(size(S,1))
    error(['DIM argument does not match channel output dimension']);
  end

  from = 'linop';
  E = S;
end



%--- linop ---
if strcmp(from,'linop') && ~strcmp(to,'linop')
  if strcmp(to,'func')
    % convert from linop to func by building a function that calls 'applychan'
    S = @(rho) applychan(S,rho,'linop',[dimA,dimB]);

  else
    % if converting to anything else, convert to choi by "gamma" involution and
    % leave rest till later (note that because of the order in which matlab
    % takes elements when reshaping, the inner reshape gives the multi-indices
    % in reverse order, i.e. E_{j,i,l,k} *not* E_{i,j,k,l}, which needs to be
    % converted to E_{jl,ik} for a choi matrix whose first subsystem is the
    % input, but because of the outer reshape we have to permute to
    % E_{l,j,k,i})
    S = reshape(permute(reshape(E,[dimB,dimB,dimA,dimA]),[1,3,2,4]),...
                dimA*dimB,dimA*dimB);
    from = 'choi';
    E = S;
  end
end


%--- choi ---
if strcmp(from,'choi')
  % convert choi to standard Choi-Jamiolkowski form, i.e. ensure reduced
  % state on first system is identity
  s = kron(TrX(E,2,[dimA,dimB])^(-1/2),eye(dimB));
  S = s*E*s;
  S = dimA*S/trace(S);
  E = S;

  if strcmp(to,'linop')
    % convert from choi to linop by "gamma" involution (E_{j,i,l,k} from inner
    % reshape becomes E_{i,k,j,l}, which the outer reshape converts to
    % E_{ki,lj})
    S = reshape(permute(reshape(E,[dimB,dimA,dimB,dimA]),[1,3,2,4]),...
                dimB^2,dimA^2);

  elseif strcmp(to,'func')
    % convert from choi to func by building a function that calls 'applychan'
    S = @(rho) applychan(S,rho,'choi',[dimA,dimB]);

  elseif ~strcmp(to,'choi')
    % convert from choi to state by taking purification; if converting to
    % anything else, convert to state and leave rest till later
    [S,dimE] = purification(E);
    from = 'state';
    E = S;
  end
end



%--- state ---
% FIXME: should convert state to "standard" state a la choi to choi
if strcmp(from,'state') && ~strcmp(to,'state')
  if strcmp(to,'func')
    % convert isometry to func by building function that calls 'applychan'
    S = @(rho) applychan(S,rho,'state',[dimA,dimB]);

  else
    % convert state to isometry by finding Schmidt decomposition for A|BE
    % partition, and "reversing" the ket on input system into a bra; if
    % converting to anything else, convert to isometry and leave rest till
    % later
    [U,D,V] = schmidt(E,[dimA,dimB*dimE]);
    S = V*diag(D)*U.';
    from = 'isom';
    E = S;
  end
end



%--- isom ---
if strcmp(from,'isom')
  if strcmp(to,'kraus')
    % convert isometry to kraus by left-multiplying by computational basis
    % state environment bras
    S = cell(dimE,1);
    for j = 1:dimE
      S{j} = E([j:dimE:dimB*dimE],:);
    end

  elseif strcmp(to,'choi') || strcmp(to,'linop')
    % convert isometry to choi by applying isometry to second half of a
    % maximally entangled state and tracing over environment
    w = zeros(dimA^2,1);
    w([1:dimA+1:dimA^2]) = 1;
    w = w*w';
    S = TrX(kron(eye(dimA),E)*w*kron(eye(dimA),E'),2,[dimA*dimB,dimE]);

    % convert choi to linop by "gamma" involution
    if strcmp(to,'linop')
      S = reshape(permute(reshape(S,[dimB,dimA,dimB,dimA]),[1,3,2,4]),...
                  [dimB^2,dimA^2]);
    end

  elseif strcmp(to,'state')
    % convert isometry to state by "reversing" the bra on the input
    % system into a ket
    S = reshape(E,dimA*dimB*dimE,1);

  elseif strcmp(to,'func')
    % convert isometry to func by building function that calls 'applychan'
    S = @(rho) applychan(S,rho,'isom',[dimA,dimB]);
  end
end


% output results
varargout{1} = S;
if nargout > 1;
   varargout{2} = [dimA,dimB,dimE];
end
