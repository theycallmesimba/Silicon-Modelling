function varargout = schmidt(psi,dim)
% SCHMIDT  Schmidt decomposition
% requires: nothing
% author: Toby Cubitt
% license: GPL2
%
%    S = SCHMIDT(PSI,DIM) returns the Schmidt coefficients of a
%    bipartite pure state PSI with subsystem dimensions specified
%    by vector DIM = [dim1,dim2].
%
%    [V1,S,V2] = SCHMIDT(PSI,DIM) in addition returns the Schmidt basis
%    vectors for the two subsystems in the columns of V1 and V2:
%
%      psi = sum_i S(i)*kron(V1(:,i),V2(:,i))


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


if nargout <= 1
  varargout{1} = svd(reshape(psi,dim(2),dim(1)).');
else
  [V1,S,V2] = svd(reshape(psi,dim(2),dim(1)).','econ');
  varargout{1} = V1;
  varargout{2} = diag(S);
  varargout{3} = conj(V2);
end
