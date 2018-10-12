% Put some useful variables in the global namespace


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


q0 = [1;0];
q1 = [0;1];
qp = (q0 + q1)/sqrt(2);
qm = (q0 - q1)/sqrt(2);
q00 = kron(q0,q0);
q01 = kron(q0,q1);
q10 = kron(q1,q0);
q11 = kron(q1,q1);
q000 = kron(q00,q0);
q001 = kron(q00,q1);
q010 = kron(q01,q0);
q011 = kron(q01,q1);
q100 = kron(q10,q0);
q101 = kron(q10,q1);
q110 = kron(q11,q0);
q111 = kron(q11,q1);

ebit = (q00 + q11)/sqrt(2);
phip = (q00 + q11)/sqrt(2);
phim = (q00 - q11)/sqrt(2);
psip = (q01 + q10)/sqrt(2);
psim = (q01 - q10)/sqrt(2);

id = eye(2);
sx = [0,1;1,0];
sy = [0,-i;i,0];
sz = [1,0;0,-1];
sp = [0,0;1,0];
sm = [0,1;0,0];

S = kron(sy,sy);

s = cell(4,1);
s{1} = id;
s{2} = sx;
s{3} = sy;
s{4} = sz;
