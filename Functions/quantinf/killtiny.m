function A = killtiny(B,varargin)

% KILLTINY   remove small values
% requires: nothing
% author: Toby Cubitt
% license: GPL2
%
%    B = killtiny(A) returns the matrix B with all real and imaginary
%    parts that are within within 1E-10 of zero set to exactly zero.
%
%    B = killtiny(A,tol) uses the threshold to TOL.


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

A = B;
idx = (abs(imag(B)) < tol);
A(idx) = real(A(idx));
idx = (abs(real(B)) < tol);
A(idx) = i*imag(A(idx));
