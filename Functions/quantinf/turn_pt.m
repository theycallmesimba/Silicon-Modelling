function [interest,x0,x1,steps] = turn_pt(last_few)


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


if (last_few(3,2)-last_few(2,2))*(last_few(2,2)-last_few(1,2)) < 0
  interest = true;
  x0 = last_few(1,1);
  x1 = last_few(3,1);
  steps = -1;
else
  interest = false;
  x0 = 0;
  x1 = 0;
  steps = 0;
end
