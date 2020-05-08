function M = tensor(varargin)

% TENSOR  Tensor product
% author: Toby Cubitt
% requires: none
% license: GPL2
%
%   m = TENSOR(a,b,c,...) returns the kronecker product of its arguments.
%
%   Each argument should either be a matrix, or a cell array containing a
%   matrix and an integer. In the latter case, the integer specifies the
%   repeat count for the matrix, e.g. TENSOR(a,{b,3},c) = TENSOR(a,b,b,b,c).


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


M = 1;
for j = 1:nargin
  if iscell(varargin{j})
    for k = 1:varargin{j}{2}
      M = kron(M,varargin{j}{1});
    end
  else
    M = kron(M,varargin{j});
  end
end
