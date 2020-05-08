function x = TrX23(p,sys,dim);

% TRX23   Partial trace of bi/tri-partite systems
% requires: nothing
% author: Toby Cubitt
% license: GPL2
%
%    X = TrX23(P,SYS,DIM) traces out system SYS of state P (a state
%    vector or a density matrix) with subsystem dimensions specified
%    by DIM.
%
%    If only one dimension is specified, i.e. DIM=[dim1], a
%    dim1 x length(p)/dim1 system is assumed.
%
%    If two are specified, i.e. DIM=[dim1,dim2], a dim1 x dim2
%    system is assumed.
%
%    DIM=[dim1,dim2,dim3] specifies a dim1 x dim2 x dim3 system
%    (duh!)


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


% check arguments
if sys > 2 && length(dim) < 3
  error('SYS greater than number of subsystems')
end
if (length(dim) == 1 && mod(length(p)/dim,1) ~= 0)...
  || length(p) ~= prod(dim)
  error('Size of P inconsistent with DIM');
end


% sort out sys and dim arguments
switch length(dim)
%  case 0
%   dim1 = 2;
%   dim2 = 1;
%   dim3 = 2;
%   if (sys == 2) sys = 3; end
 case 1
  dim1 = dim(1);
  dim2 = 1;
  dim3 = length(p)/dim1;
  if (sys == 2) sys = 3; end
 case 2
  dim1 = dim(1);
  dim2 = 1;
  dim3 = dim(2);
  if (sys == 2) sys = 3; end
 case 3
  dim1 = dim(1);
  dim2 = dim(2);
  dim3 = dim(3);
end


% calculate partial trace
switch any(size(p)==1)
 case 1
  % state vector
  if size(p,1) == 1
    p = p';
  end

  switch sys
   case 1
    x = reshape(p,dim2*dim3,dim1);
   case 2
    x = reshape(permute(reshape(p,dim3,dim2,dim1),[1,3,2]),dim1*dim3,dim2);
   case 3
    x = reshape(p,dim3,dim1*dim2).';
  end
  x = x*x';


 case 0
  % density matrix
  switch sys
   case 1
    x=zeros(dim2*dim3);
    indx=(1:dim2*dim3);
    for k=0:dim1-1
      x=x+p(indx+k*dim2*dim3,indx+k*dim2*dim3);
    end

   case 2
    x=zeros(dim1*dim3);
    indx=kron(ones(1,dim1),[1:dim3]) + ...
	 kron(dim2*dim3*[0:dim1-1],ones(1,dim3));
    for k=0:dim2-1
      x=x+p(indx+k*dim3,indx+k*dim3);
    end

   case 3
    x=zeros(dim1*dim2);
    indx=dim3*(0:dim1*dim2-1);
    for k=1:dim3
      x=x+p(indx+k,indx+k);
    end
  end
end
