function [x, obj, Y, info] = runLMI(c, F, G, varargin)

% RUNLMI  Solve a semi-definite program given in dual (linear matrix
%         inequality) form
% requires: SeDuMi
% author: Toby Cubitt
% license: GPL2
%
%     [x, OBJ, Y, INFO] = runLMI(c, F, G) solves the following SDP
%     (where >= denotes matrix inequality):
%
%       minimize     c.'*x
%       subject to   sum_i F{k,i}*x(i) >= G{k}
%
%     [x, OBJ, Y, INFO] = runLMI(c, F, G, A, b) solves the following SDP
%     with additional equality constraints:
%
%       minimize     c.'*x
%       subject to   sum_i F{k,i}*x(i) >= G{k}
%                    A*x == b
%
%     The vector c defines the linear objective function, F is a
%     2-dimensional cell array containing the LMI coefficient matrices, G is
%     a 1-dimensional cell array containing the LMI constant matrices,
%     matrix A and vector b specify the equality constraints.
%
%     x is the optimal solution vector, OBJ is the value of the objective
%     function at that point, Y is the dual solution matrix. INFO is a
%     SeDuMi structure containing information about the outcome of the
%     numerical optimisation.
%
%     [x, OBJ, Y, INFO] = runLMI(c, F, G, 'real') or
%     [x, OBJ, Y, INFO] = runLMI(c, F, G, A, b, 'real') assumes F and G are
%     real symmetric matrices. Otherwise, G and F can be Hermitian.
%
%     If your problem only requires a single LMI, G is allowed to be a
%     single matrix rather than a cell array (F will be a
%     1-dimensional cell array in this case).
%
%
%     See also RUNSDP, SEDUMI.
%
%     Note that the optimisation problems solved by RUNLMI and RUNSDP are
%     not strictly dual to each other: if [x,f,Y] = runLMI(c,F,G), then
%     [Y,-f,-x] = runSDP(-G,F,c).


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


% sort out input arguments
if isnumeric(G)
  if size(F,2) == 1
    F = F';
  end
  if size(F,1) ~= 1
    error('inconsistent numbers of matrices in F and G');
  end
  G = {G};
end

R = 'complex';
if length(varargin) >= 3
  R = varargin{3};
elseif length(varargin) >= 1 && ischar(varargin{1})
  R = varargin{1};
end


% get dimensions etc.
nvars = length(c);
nLMI = size(F,1);
dim = zeros(nLMI,1);
for k = 1:nLMI
  dim(k) = size(F{k,1},1);
end



%--- equality constraints ---
At = sparse([]);
bt = sparse([]);
if length(varargin) >= 2
  At = sparse(varargin{1});
  bt = sparse(varargin{2});
  K.f = size(At,1);
end


%--- LMI constraints ---
K.s = dim;
if ~ (strcmp(R,'real') || strcmp(R,'re'))
  K.scomplex = [1:nLMI];
end
for k = 1:nLMI
  Atmp = sparse([]);
  for j = 1:nvars
    % !!! -1 because SeDuMi has LMI in form: G - x*F >= 0 !!!
    Atmp = [Atmp,vec(-sparse(F{k,j}))];
  end
  At = [At;Atmp];
  % !!! -1 because SeDuMi has LMI in form: G - x*F >= 0 !!!
  bt = [bt;vec(-sparse(G{k}))];
end


%--- objective function ---
% !!! -1 because SeDuMi maximizes dual rather than minimizing !!!
c = -sparse(c);
% ensure c is a column vector
if size(c,1) == 1
  c = c.';
end


%% %--- display problem definition ---
%% disp('objective function: c');
%% disp(full(c));
%% disp('');

%% A = [];
%% if length(varargin) >= 2
%%    A = varargin{1};
%% end

%% disp('LMI:')
%% disp('F');
%% for j = 1:nvars
%%   disp(full(mat(At(size(A,1)+1:size(A,1)+dim(1)^2,j))));
%% end
%% for k = 2:nLMI
%%   for j = 1:nvars
%%     disp('');
%%     disp(full(mat(At(size(A,1)+sum(dim(1:k-1).^2)+1:...
%% 		     size(A,1)+sum(dim(1:k).^2),j))));
%%   end
%% end
%% disp('');

%% disp('G');
%% disp(full(mat(bt(size(A,1)+1:size(A,1)+dim(1)^2))));
%% for k = 2:nLMI
%%   disp('');
%%   disp(full(mat(bt(sum(size(A,1)+dim(1:k-1).^2)+1:...
%% 		  size(A,1)+sum(dim(1:k).^2)))));
%% end
%% disp('');

%% if ~ (A == [])
%%   disp('equality constraint:')
%%   disp('A')
%%   disp(full(At(1:size(A,1),:)));

%%   disp('b');
%%   disp(full(bt(1:size(A,1))));
%% end


%--- run SDP ---
[x,y,info] = sedumi(At, c, bt, K);
Y{1} = mat(x(1:dim(1)^2));
for k = 2:nLMI
  Y{k} = mat(x(sum(dim(1:k-1).^2)+1:sum(dim(1:k).^2)));
end
x = y;
obj = full(-c.'*x);  % undo sign-change introduced into objective function
