function P = Psym(d)

% PSYM  Projector onto symmetric subspace
% requires: nothing
% author: Toby Cubitt
% License: GLP2
%
%   Psym(d) returns the projector onto the symmetric subspace of a d x d
%   system.

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


P = zeros(d^2);

for j = 1:d
  a = zeros(d);
  a(j,j) = 1;
  P = P + kron(a,a);
end

for j = 1:d
  for k = j+1:d
    a = zeros(d,1); a(j) = 1;
    b = zeros(d,1); b(k) = 1;
    a = (kron(a,b) + kron(b,a));
    a = a*a'/2;
    P = P + a;
  end
end
