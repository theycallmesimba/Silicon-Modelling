function varargout = maxEfrac(rho)
% maxEfrac   Maximum entangled fraction
% requires: unpauli.m
% author: Toby Cubitt
% license: GPL2
%
%    S = maxEfrac(RHO) returns the maximum entangled fraction (also
%    known as maximum singlet fraction) of 2x2 density matrix RHO,
%    defined as the
%
%      max phi'*rho*phi
%
%    where phi is a two-qubit maximally entangled state
%
%    [O1,S,O2] = maxEfrac(RHO) in addition returns the orthogonal
%    matrices that achieve the maximum in the Hilbert-Schmidt
%    representation.


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


R = pauli(rho);
T = abs(R(2:4,2:4)).*sign(real(R(2:4,2:4)));

[O1,S,O2] = svd(T);

S = diag([1,1,det(O1)])*S*diag([1,1,det(O2)]);
O1 = O1*diag([1,1,det(O1)]);
O2 = O2*diag([1,1,det(O2)]);

Q = (O2*diag([1,1,-1])*O1.').';
P(1,1) = 1;
P(2:4,2:4) = Q;

phi = unpauli(P);

varargout{1} = trace(phi*rho);
if nargout >= 2
  varargout{2} = phi;
end
if nargout == 4
  varargout{3} = O1;
  varargout{4} = O2;
end
