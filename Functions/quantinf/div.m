function x = div(a,b)

% DIV  Integer division
% requires: nothing
% author: Toby Cubitt
% license: GPL2
%
%     X=DIV(A,B) returns the integer part of a/b,
%
%                     {  floor(a./b)   a*b >= 0
%     i.e. div(a,b) = {
%                     {  ceil(a./b)    a*b <= 0


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


% scalar a and b
if (length(a) == 1) && (length(b) == 1)
  if a*b >= 0
    x = floor(a./b);
  else
    x = ceil(a./b);
  end
else
  pos = find(a.*b >= 0);
  neg = find(a.*b < 0);
  if length(a) == 1
    x = zeros(size(b));
    x(pos) = floor(a./b(pos));
    x(neg) = ceil(a./b(neg));
  elseif length(b) == 1
    x = zeros(size(a));
    x(pos) = floor(a(pos)./b);
    x(neg) = ceil(a(neg)./b);
  else
    x = zeros(size(a));
    x(pos) = floor(a(pos)./b(pos));
    x(neg) = ceil(a(neg)./b(neg));
  end
end
