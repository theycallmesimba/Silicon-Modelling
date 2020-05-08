function pprintm(M,varargin)

% PPRINTM   Matrix pretty-printer
% requires: nothing
% author: Toby Cubitt
% license: GPL2
%
%    pprintm(M) pretty-prints a matrix M in a format that makes it easier to
%    see the structure of the matrix than the default output format does.
%
%    pprintm(M,LABELS) pretty-prints the matrix, using LABELS to label the
%    rows and columns (see below).
%
%    pprintm(M,ROWLABELS,COLLABELS) pretty-prints the matrix, using ROWLABELS
%    and COLLABELS to label the rows and columns, respectively (see below).
%
%    If M is a row or column vector rather than a matrix, only column or
%    row labels (respectively) will be displayed.
%
%    When specified, LABELS, ROWLABELS and COLLABELS can take the following
%    forms:
%
%      column vector of strings - strings used as labels, cyclically
%      cell array {d,n}         - labels for n qudits
%      vector of integers       - labels for qudits with specified dims


%% Copyright (C) 2011-2012 Toby Cubitt
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



% generate any pre-defined labels
for l = 1:length(varargin)
  if length(varargin) >= l
    % labels for n qudits
    if iscell(varargin{l})
      d = varargin{l}{1};
      n = varargin{l}{2};
      for j = 0:d^n-1
	lbl = "|";
	for k = n-1:-1:0
	  lbl = strcat(lbl, sprintf("%d",mod(floor(j/d^k),d)));
	end
	lbl = strcat(lbl,">");
	labels{l}(j+1,:) = lbl;
      end

    % labels for qudits with specified dimensions
    elseif isnumeric(varargin{l}(1))
      d = varargin{l};
      n = length(d);
      for j = 0:prod(d)-1
	lbl = "|";
	for k = 1:n
	  lbl = strcat(lbl, sprintf("%d",mod(floor(j/prod(d(k+1:end))),d(k))));
	end
	lbl = strcat(lbl,">");
	labels{l}(j+1,:) = lbl;
      end

    else
	labels{l} = varargin{l};
    end

    varargin{l} = labels{l};
  end
end



% get length of longest matrix entry
maxlen = 0;
for j=1:size(M)(1)
  for k=1:size(M)(2)
    if isreal(M(j,k))
      maxlen = max(maxlen,length(sprintf("%.3g",M(j,k))));
    elseif isreal(i*M(j,k))
      maxlen = max(maxlen,length(sprintf("%.3gi",M(j,k))));
    else
      maxlen = max(maxlen,length(sprintf("%.3g + %.3gi", ...
					 real(M(j,k)),imag(M(j,k)))));
    end
  end
end

% column labels
fprintf("\n")
c = length(varargin);
if c > 0 && size(M)(2) > 1
  colfmt = sprintf("%%%ds  ", maxlen);
  % (note that all strings in a vector of strings have the same length, so
  % length of all row labels is the same)
  rlen = length(varargin{1}(1,:));
  if size(M)(1) > 1; fprintf(blanks(rlen+3)); end
  for j=1:size(M)(2)
    fprintf(colfmt, varargin{c}(mod(j-1,size(varargin{c})(1))+1,:));
  end
  fprintf("\n\n");
  % if column labels are longer than longest matrix entry, pad matrix entries
  maxlen = max(maxlen,length(varargin{c}(1,:)));
end


% matrix entry format
fmt = sprintf("%%%ds  ",maxlen);

for j=1:size(M)(1)
  % row labels
  if length(varargin) > 0 && size(M)(1) > 1
    fprintf("%s   ", varargin{1}(mod(j-1,size(varargin{1})(1))+1,:));
  end
  % matrix rows
  for k=1:size(M)(2)
    if isreal(M(j,k))
      fprintf(fmt,sprintf("%.3g",M(j,k)));
    elseif isreal(i*M(j,k))
      fprintf(fmt,sprintf("%.3gi",imag(M(j,k))));
    else
      if imag(M(j,k)) < 0
	if imag(M(j,k)) == -1
	  fprintf(fmt,sprintf("%.3g - i",real(M(j,k))));
	else
	  fprintf(fmt,sprintf("%.3g - %.3gi",real(M(j,k)),-imag(M(j,k))));
	end
      else
	if imag(M(j,k)) == 1
	  fprintf(fmt,sprintf("%.3g + i",real(M(j,k))));
	else
	  fprintf(fmt,sprintf("%.3g + %.3gi",real(M(j,k)),imag(M(j,k))));
	end
      end
    end
  end
  fprintf("\n")
end

fprintf("\n")
