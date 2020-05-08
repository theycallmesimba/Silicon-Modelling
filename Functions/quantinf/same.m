function b = same(a,b,varargin)

% SAME   test whether two matrices are identical
% requires: nothing
% author: Toby Cubitt
% license: GPL2
%
%    B = SAME(A,B) tests whether matrices A and B are identical to within
%    an entry-wise tolerance of 1E-10.
%
%    B = SAME(A,B,TOL) tests whether matrices A and B are identical to
%    within an entry-wise tolerance of TOL.


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


if isempty(varargin)
  tol = 1E-10;
else
  tol = varargin{1};
end

if exist('full')
  b = full(all(all(abs(a - b) < tol)));
else
  % octave
  b = all(all(abs(a - b) < tol));
end
