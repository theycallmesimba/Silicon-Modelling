function [X, obj, y, info] = runSDP(C, A, b, varargin)

% RUNSDP  Solve a semi-definite program given in primal form
% requires: SeDuMi
% author: Toby Cubitt
% license: GPL2
%
%     [X, OBJ, y] = runSDP(C, A, b) solves the following SDP
%     (where >= denotes matrix inequality):
%
%       minimize     trace(C*X)
%       subject to   X >= 0
%                    trace(A{i}*X) == b(i)
%
%     The matrix C defines the objective function, A is a 1-dimensional cell
%     array containing the matrices appearing in the equality constraints, b
%     is a vector containing the corresponding constants.
%
%     X is the optimal solution matrix, OBJ is the value of the objective
%     function at that point, y is the dual solution vector. INFO is a
%     SeDuMi structure containing information about the outcome of the
%     numerical optimisation.
%
%     [X, OBJ, y] = runSDP(C, A, b, 'real') restricts X to be a real
%     symmetric matrix. Otherwise, X can be Hermitian.
%
%     If your problem only requires a single constraint, A is allowed to be
%     a single matrix rather than a cell array (b will be a scalar in this
%     case).
%
%
%     See also RUNLMI, SEDUMI.
%
%     Note that the optimisation problems solved by RUNSDP and RUNLMI are
%     not strictly dual to each other: if [X,f,y] = runSDP(C,A,b), then
%     [-y,-f,X] = runLMI(b,A,-C).


%% Copyright (C) 2011 Toby Cubitt
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


% sort out input arguments
if isnumeric(A)
  A = {A};
end
if ~ length(b) == length(A);
  error('dimensions of A and b are inconsistent');
end

R = 'hermitian';
if length(varargin) >= 1
  R = varargin{1};
end


% get dimensions etc.
dim = size(C,1);
K.s = dim;
if ~ (strcmp(R,'real') || strcmp(R,'re'))
  K.scomplex = [1];
end


%--- equality constraints ---
At = sparse([]);
for k=1:length(A)
    At = [At;vec(sparse(A{k})).'];
end


%--- objective function ---
C = vec(sparse(C));


%--- run SDP ---
[X,y,info] = sedumi(At, b, C, K);
X = mat(X);
obj = full(trace(mat(C)*X));
