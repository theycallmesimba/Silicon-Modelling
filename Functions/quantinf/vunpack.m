function varargout = vunpack(v,varargin)

% VUNPACK   unpack vector into variables
% requires: nothing
% author: Toby Cubitt
% license: GPL2
%
%    [s,...,v,...,A,...] = vunpack(V,s,...,v,...,A,...)
%    unpacks vector V into variables of types specified
%    by the input arguments. These can be a mixture of
%    scalars, vectors and arrays. Elements of arrays are
%    assigned in row-major order. Note that the contents
%    of the specification input variables (all input
%    arguments except V) are ignored.


%% Copyright (C) 2004-2011 Toby Cubitt
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


for j = 1:nargin-1
  s = prod(size(varargin{j}));
  varargout{j} = reshape(v(1:s),size(varargin{j}));
  v = v(s+1:end);
end
