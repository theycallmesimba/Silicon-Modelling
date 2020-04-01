function x = Tx23(p,sys,dim)

% Tx   Partial transpose of bi/tri-partite systems
% requires: nothing
% author: Toby Cubitt
% license: GPL2
%
%    X = Tx23(RHO,SYS,DIM) gives the partial transpose of a matrix RHO
%    with respect to subsystem SYS, where the subsystem dimensions are
%    specified by vector DIM
%
%    If no dimensions are specified, ie DIM = [], Tx assumes a 2x2
%    system.
%
%    If only one is specified, i.e. DIM=[dim1], it assumes a
%    dim1 x dim1 system.
%
%    If two dimensions are specified, i.e. DIM=[dim1,dim2],
%    dim1 x dim2 is assumed.
%
%    DIM=[dim1,dim2,dim3] specifies three subsystems (duh!)


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


switch length(dim)
 case 0
  dim1 = 2;
  dim2 = 1;
  dim3 = 2;
 case 1
  dim1 = dim(1);
  dim2 = 1;
  dim3 = dim1;
 case 2
  dim1 = dim(1);
  dim2 = 1;
  dim3 = dim(2);
 case 3
  dim1 = dim(1);
  dim2 = dim(2);
  dim3 = dim(3);
end
if (dim2 == 1 & sys == 2) sys = 3; end


switch sys
 case 1
  x = p;
  indx = (1:dim2*dim3);
  for j=0:dim1-1
    for k = 0:dim1-1
      if j==k continue; end
      x(indx+j*dim2*dim3,indx+k*dim2*dim3) = ...
	  p(indx+k*dim2*dim3,indx+j*dim2*dim3);
    end
  end

 case 2
  x = p;
  indx=kron(ones(1,dim1),[1:dim3]) + ...
       kron(dim2*dim3*[0:dim1-1],ones(1,dim3));
  for j = 0:dim2-1
    for k = 0:dim2-1
      if j==k continue; end
      x(indx+j*dim3,indx+k*dim3) = ...
	  p(indx+k*dim3,indx+j*dim3);
    end
  end

 case 3
  x = p;
  indx = (1:dim3);
  for j=0:dim1*dim2-1
    for k = 0:dim1*dim2-1
      x(indx+j*dim3,indx+k*dim3) = ...
          p(indx+j*dim3,indx+k*dim3)';
    end
  end

end
