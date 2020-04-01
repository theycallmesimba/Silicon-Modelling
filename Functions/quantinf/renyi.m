function S = renyi(p,rho)

% RENYI   Renyi entropy
% requires: nothing
% author: Toby Cubitt
% license: GPL2
%
%    S = RENYI(P,RHO) returns the Renyi p-entropy of state RHO (a
%    probability distribution or a density matrix)


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


if size(rho,1) == 1 || size(rho,2) == 1
  if size(rho,1) == 1
    rho = rho.';
  end
  switch p
   case 1,
    S = -rho.'*log2(rho+(rho==0));
   case 0,
    tol = length(rho)*eps(norm(rho));
    S = log2(sum(rho>tol));
   otherwise
    S = 1/(1-p) * log2(sum(rho.^p));
  end

else
  switch p
   case 1,
    e = eig(rho);
    S = -e'*log2(e+(e==0));
   case 0,
    S = log2(rank(rho));
   otherwise,
    S = 1/(1-p) * log2(trace(rho^p));
  end
end
