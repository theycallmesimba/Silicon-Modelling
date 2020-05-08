function v = vpack(varargin)

% VPACK   Pack variables into a vector
% requires: nothing
% author: Toby Cubitt
% license: GPL2
%
%    vpack(s,...,v,...,A,...) packs its arguments into a
%    single output vector. Arguments can be any mixture of
%    scalars, vectors, and arrays. Arguments are packed in
%    order, with elements of arrays packed in row-major order.


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


v = [];
for j = 1:length(varargin)
  v = [v,varargin{j}(1:prod(size(varargin{j})))];
end
