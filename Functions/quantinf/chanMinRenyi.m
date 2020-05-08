function Smin = chanMinRenyi(p,E,rep,dim,varargin)

% CHANMINRENYI  Minimum output Renyi entropy
% requires: renyi, applychan, randPsi
% author: Toby Cubitt
% license: GPL2
%
%    Smin = chanMinRenyi(p,E,REP,DIM) numerically minimizes the output
%    Renyi p-entropy of channel E, whose representation is specified
%    by REP and whose input and output dimensions dimA and dimB are
%    given by the vector DIM = [dimA dimB].
%
%    Smin = minRenyi(p,E,REP,DIM,PSI0) uses PSI0 as the initial point
%    for the optimization.
%
%    Smin = minRenyi(p,E,REP,DIM,PSI0,TRIES) takes the best result
%    over TRIES runs of the numerical optimization.
%
%    Smin = minRenyi(p,E,REP,DIM,PSI0,TRIES,OPTS) passes OPTS as
%    options to the optimization routine (see OPTIMSET).
%
%    To leave one of the optional arguments unspecified, set it to [].


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


%--- sort out input arguments ---
tries = 1;
opts = optimset();
if length(varargin) >= 1 && ~isempty(varargin{1})
  psi0 = varargin{1};
  psi0 = [real(psi0),imag(psi0)];
end
if length(varargin) >= 2 && ~isempty(varargin{2})
  tries = varargin{2};
end
if length(varargin) >= 3 && ~isempty(varargin{3})
  opts = varargin{3};
end


%--- minimization functions ---
fcon = @(psi) renyi(p, applychan(E,psi(:,1)+i*psi(:,2),rep,dim));
con = @(psi) deal(0, norm(psi(:,1)+i*psi(:,2)) - 1);
func = @(psi) renyi(p, applychan(E,normalise(psi(:,1)+i*psi(:,2)),rep,dim));


%--- do minimization ---
% do first optimisation run
flag = -2;
j = 0;
while flag < 0 && j < tries
  try
    if length(varargin) < 1 || isempty(varargin{1})
      psi0 = randPsi(dim(1));
      psi0 = [real(psi0),imag(psi0)];
    end
    [psi1,r,flag] = fminunc(func,psi0,opts);
    %[psi1,r,flag] = fmincon(fcon,psi0,[],[],[],[],[],[],con);
    r = real(r);
    if flag < 0
      j = j + 1;
    end
  catch
    j = j + 1
  end
end


% if entropy is smaller at explicitly supplied starting point, use that
if length(varargin) >= 1 && ~isempty(varargin{1})
  Smin = real(renyi(p,applychan(E,psi0(:,1)+i*psi0(:,2),rep,dim)));
  if r < Smin
    Smin = r;
  end
else
  Smin = r;
end


% do rest of optimization runs
for k = 2:tries
  flag = -2;
  j = 0;
  while flag < 0 && j < tries
    try
      psi0 = randPsi(dim(1));
      psi0 = [real(psi0),imag(psi0)];
      [psi1,r,flag] = fminunc(func,psi0,opts);
      %[psi1,r] = fmincon(fcon,psi0,[],[],[],[],[],[],con);
      r = real(r);
      if r < Smin
        Smin = r;
      end
      if flag < 0
        j = j + 1;
      end
    catch
      j = j + 1;
    end
  end
end
