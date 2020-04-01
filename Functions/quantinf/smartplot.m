function varargout = ...
    smartplot(plot_func,x0,x1,steps,interest_func,last_n,varargin)
% SMARTPLOT   Generate samples of a scalar function for plotting
% requires: nothing
% author: Toby Cubitt
% license: GPL2
%
%    [abscissa,ordinate,[exitflag]] =
%       SMARTPLOT(plot_func,x0,x1,steps,
%                 interest_func,last_n,[options,[P1,...]])
%
%    SMARTPLOT returns two vectors [abscissa,ordinate] containing x,y
%    pairs, with x scanned over the interval x0 to x1, for the
%    function supplied in PLOT_FUNC. These can be passed directly to
%    Matlab's PLOT function to plot a graph of PLOT_FUNC.
%
%    SMARTPLOT optionally also returns an exit flag. If EXITFLAG is:
%      0   then SMARTPLOT successfully calculated all points
%      > 0 then SMARTPLOT had to skip this many points because
%          PLOT_FUNC indicated errors via it's ERROR output
%      < 0 then SMARTPLOT encountered a fatal error and tried to
%          output all the points it did manage to calculate, of
%          which there should be this many
%
%
%    SMARTPLOT initially starts sampling PLOT_FUNC from x0, in
%    increments such that the interval x0 to x1 is covered by STEPS
%    samples. If the user-supplied function INTEREST_FUNC indicates
%    that the interval spanned by the last few points contains an
%    interesting feature (e.g. maximum, zero, etc.), that interval is
%    sampled at a higher resolution (i.e. in smaller increments)
%    before continuing the sampling between x0 and x1.
%
%    The procedure is recursive; if an interesting feature is detected
%    in an interval being scanned at a higher resolution, that
%    interval is sampled at an even higher resolution. The recursion
%    is limited to some maximum 'depth', set in the options (see
%    OPTIONS below).
%
%
%    PLOT_FUNC must take at least one real, scalar input x. It may
%    also take any number of further fixed input parameters, whose
%    values must then be supplied in optional arguments P1,...  It
%    should return the scalar value y corresponding to x. It may
%    optionally return a second boolean argument ERRORFLAG which is
%    set to true if it could not calculate the current
%    point, causing SMARTPLOT to skip it.
%
%
%    X0,X1 define the interval over which x should be varied. Points
%    at the ends of the interval are included. Thus if no interesting
%    features are detected, SMARTPLOT will calculate STEPS points
%    spaced evenly over the interval X0,X1 with the first and last
%    points at X0 and X1 respectively. The spacing between points will
%    therefore be (X1-X0)/(STEPS-1).
%
%
%    INTEREST_FUNC must take one n x 2 array as input, containing the
%    x,y values of the last few points sampled. The number of points
%    to pass is given by the LAST_N argument to
%    SMARTPLOT. INTEREST_FUNC must return four values:
%
%      [interest (boolan), start (scalar), end (scalar), steps (int)]
%
%      INTEREST indicates whether or not there is an interesting
%      feature in the interval spanned by the points passed to
%      INTEREST_FUNC.
%
%      START returns the start of the interval to scan at a higher
%      resolution. If START is set to the special string 'start', the
%      start of the interval is set to the value of the X0 argument to
%      SMARTPLOT.
%
%      END returns the end of the interval. Usually the interval
%      should contain the feature of interest. If END is set to the
%      special string 'end', the end of the interval is set to the
%      value of the X1 argument to SMARTPLOT.
%
%
%      STEPS returns the number of points to sample from the interval
%      BEGIN to END. If STEPS is 0, the default set by
%      OPTIONS.REFINE_STEPS will be used (see below). If STEPS is
%      negative, SMARTPLOT will calculate the number of steps
%      required to put (ABS(STEPS) * OPTIONS.REFINE_STEPS) points
%      between each point at the previous resolution.
%
%
%    Optional arguement OPTIONS, if supplied, must contain a structure
%    over-riding the default parameters. The parameters are:
%
%      OPTIONS.REFINE_DEPTH sets the maximum depth of the recursive
%      refinement of intervals containing interesting
%      features. I.e. when an interesting feature is detected,
%      SMARTPLOT will only increase the sample resolution this many
%      times. Note: setting OPTIONS.REFINE_DEPTH to 0 will allow
%      unlimited recursion, which could lead to an infinite loop.
%
%      OPTIONS.REFINE_STEPS sets the default number of steps to
%      sample from the interval returned by INTEREST_FUNC. This can
%      be overridden by the STEPS output from INTEREST_FUNC.
%
%      OPTIONS.EXPENSIVE_FUNCTION tells SMARTPLOT whether PLOT_FUNC is
%      expensive time-wise to evaluate. If set to true, SMARTPLOT will
%      save points it would otherwise throw away, and try to avoid
%      calculating duplicate points (see DUPLICATE_THRESHOLD
%      below). It is recommended this be set to TRUE.
%
%      OPTIONS.DUPLICATE_THRESHOLD sets the threshold distance between
%      points below which SMARTPLOT considers them the same
%      point. Only relevant if EXPENSIVE_FUNCTION is set to true.
%      Default is the interval that there would be between points at a
%      recursion depth 1 more than the default maximum, assuming
%      REFINE_STEPS steps at each depth.
%
%      OPTIONS.GUESS_NUM_POINTS and OPTIONS.GUESS_NUM_SAVED set
%      estimated numbers of points that SMARTPLOT will have to
%      calculate or save (see EXPENSIVE_FUNCTION above). SMARTPLOT
%      will not fail if these numbers are exceeded, but it could have
%      a speed impact.
%
%      OPTIONS.ENABLE_FUN_CHECK sets whether to check that the
%      PLOT_FUN and INTEREST_FUN functions take and return the
%      correct number and types of arguments, before starting. If
%      SMARTPLOT is being used in script, it is recommended this be
%      set to FALSE once the script is debugged.
%
%      PLOTFUN_ERRORFLAG sets whether PLOTFUN returns the optional
%      second ERRORFLAG argument or not. Ignored if ENABLE_FUN_CHECK
%      is TRUE, since SMARTPLOT can then determine this itself.
%
%      OPTIONS.LAST_FEW_TEST sets the array of x,y values to pass to
%      INTEREST_FUN when testing it. Only used if ENABLE_FUN_CHECK is
%      TRUE.


%% Copyright (C) 2004-2010 Toby Cubitt
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



% use default options if none specified
if (length(varargin) == 0) | isempty(varargin{1})
  options = struct('refine_depth',3, ...
		   'refine_steps',10, ...
		   'expensive_function',true, ...
		   'duplicate_threshold',(x1-x0)/(steps-1)/(10^3), ...
		   'guess_num_points',10*steps, ...
		   'guess_num_saved',5, ...
		   'enable_fun_check',true, ...
		   'plotfun_errorflag',true, ...
		   'last_few_test',ones(last_n,2));
else
  options = varargin{1};
end
% any further optional arguments are plot_func parameters
func_input = cell(1);
func_input{1} = x0;
func_input{2:length(varargin)} = varargin{2:end};
if options.plotfun_errorflag
  func_nargout = 2;
else
  func_nargout = 1;
end


% clean-up input arguments
if (x0 > x1) tmp1 = x0; x0 = x1; x1 = tmp1; end


% check input arguments
if (rem(steps,1) ~= 0) | (steps <= 0)
  error('MATLAB:badopt',['Argument STEPS to SMARTPLOT must be a' ...
	' positive integer']);
end

if (rem(last_n,1) ~= 0) | (last_n < 1)
  error('MATLAB:badarg',['Argument LAST_N to SMARTPLOT must ' ...
	' be an integer greater than 0']);
end

if (rem(options.refine_steps,1) ~= 0) | options.refine_steps < 0
  error('MATLAB:badopt',['OPTIONS.refine_steps must be a' ...
	' positive integer']);
end

if (rem(options.refine_depth,1) ~= 0) | (options.refine_depth < 0)
  error('MATLAB:badopt',['OPTIONS.refine_depth must be a positive' ...
	' integer or 0']);
end

if options.enable_fun_check
  try
    [tmp1,tmp2] = feval(plot_func,func_input{:});
    func_nargout = 2;
  catch
    err = lasterror;
    if strcmp(err.identifier,'MATLAB:maxlhs') | ...
	  strcmp(err.identifier,'MATLAB:TooManyOutputs')
      try
	tmp1 = feval(plot_func,func_input{:});
	func_nargout = 1;
      catch
	disp(lasterr);
	error('MATLAB:badarg',['Argument PLOT_FUNC to SMARTPLOT must' ...
	      ' resolve to a function that takes a number of' ...
	      ' inputs equal to 1 plus the number of fixed' ...
	      ' parameters, and returns either one scalar output or' ...
	      ' two outputs: [scalar,boolean]']);
      end
    else
      disp(lasterr);
      error('MATLAB:badarg',['Argument PLOT_FUNC to SMARTPLOT must' ...
	    ' resolve to a function that takes a number of' ...
	    ' inputs equal to 1 plus the number of fixed' ...
	    ' parameters, and returns either one scalar output or' ...
	    ' two outputs: [scalar,boolean]']);
    end
  end

  try
    [tmp1,tmp2,tmp3,tmp4] = feval(interest_func,options.last_few_test);
  catch
    disp(lasterr);
    error('MATLAB:badarg',['Argument INTEREST_FUNC to SMARTPLOT' ...
	  ' must resolve to a function that takes one n x 2 input' ...
	  ' vector and returns four outputs:' ...
	  ' [boolean,scalar,scalar,integer]']);
  end
end

clear tmp*;



% initialization
abscissa = zeros(options.guess_num_points,1);
ordinate = zeros(options.guess_num_points,1);
abscissa(:) = NaN;
ordinate(:) = NaN;
last_few = zeros(last_n,2);
last_few(:) = NaN;
saved_points = zeros(options.guess_num_saved,2);
saved_points(:) = NaN;
saved_n = 0;
refine_interval = zeros(options.refine_depth,2);
refine_res = zeros(options.refine_depth,1);
refine_steps = zeros(options.refine_depth,1);
refine_indx = zeros(options.refine_depth,1);
exitflag = 0;
errorflag = 0;

indx = 1;
refine_level = 1;
refine_interval(1,:) = [x0,x1];
refine_steps(1) = steps-1;
refine_indx(1) = 0;
refine_res(1) = (x1-x0)/(steps-1);



% the function being plotted could be very expensive, so if an
% error occurs try not to lose everything
% [yes, this happened to me once before I added this <sigh> - tsc]
try

  % loop until finished (!)
  while refine_level > 0;

    % set next abscissa value to caluclate
    x = refine_interval(refine_level,1) + ...
	refine_indx(refine_level)*refine_res(refine_level);

    % check if any previously saved points can be inserted (only
    % active if expensive_function option set since if it's no, no
    % points are ever saved)
    find_indx = find((saved_points(:,1) < x - options.duplicate_threshold));
    if length(find_indx) > 0
      tmp = sort(saved_points(find_indx,:),1);
      abscissa(indx:indx+length(find_indx)-1) = squeeze(tmp(:,1));
      ordinate(indx:indx+length(find_indx)-1) = squeeze(tmp(:,2));
      saved_n = saved_n - length(find_indx);
      saved_points(find_indx,:) = [];
      indx = indx + length(find_indx);
    end


    % if we've already calculated one or more points near the current
    % point, don't waste time calculating it, instead retreive the
    % nearby points from the array of saved ones (only active if
    % expensive_function option set since if it's not, no points are
    % ever saved)
    find_indx = find(abs(saved_points(:,1) - x) <= ...
		     options.duplicate_threshold);
    if length(find_indx) > 0
      tmp = sort(saved_points(find_indx,:),1);
      abscissa(indx:indx+length(find_indx)-1) = squeeze(tmp(:,1));
      ordinate(indx:indx+length(find_indx)-1) = squeeze(tmp(:,2));
      saved_n = saved_n - length(find_indx);
      saved_points(find_indx,:) = [];
      indx = indx + length(find_indx) - 1;

    else

      % calculate next point
      abscissa(indx) = x;
      func_input{1} = x;
      if func_nargout == 2
	[ordinate(indx),errorflag] = feval(plot_func,func_input{:});
      else
	ordinate(indx) = feval(plot_func,func_input{:});
      end

      % if plot_func returned an error, skip to next point
      if errorflag
	exitflag = exitflag + 1;
	refine_indx(refine_level) = refine_indx(refine_level) + 1;
	indx = indx - 1;

      else

% 	% update vector storing last few points
% 	last_few = circshift(last_few,-1);
% 	last_few(end,:) = [abscissa(indx),ordinate(indx)];

	% if we've got enough points to pass to interest_func
	if (indx >= last_n)
	  last_few = [abscissa(indx-last_n+1:indx), ...
		      ordinate(indx-last_n+1:indx)];
	  [interesting,x0,x1,steps] = feval(interest_func,last_few);

	  % if we've been told to abort current refinement level, and
          % we're not at the lowest refinement level, drop down a
          % level
	  if (interesting < 0) & (refine_level ~= 1)
	    refine_level = refine_level - 1;


	  % if we've found an interesting point and we're not already
          % at the maximum refinement level, go up a level
	  elseif (interesting > 0) & (refine_level ~= options.refine_depth)
	    refine_level = refine_level + 1;
	    if (x0 == 'start') x0 = refine_interval(1,1); end
	    if (x1 == 'end') x1 = refine_interval(1,2); end
	    refine_interval(refine_level,:) = [x0,x1];
	    if steps == 0
	      steps = options.refine_steps + 1;
	    elseif steps < 0
	      steps = round(abs(steps) * ...
			    round((x1-x0)/refine_res(refine_level-1)) * ...
			    options.refine_steps + 1);
	    end
	    refine_steps(refine_level) = steps;
	    refine_res(refine_level) = (x1-x0)/(steps-1);
	    refine_indx(refine_level) = -1; % it gets incremented below

	    % rewind back to beginning of refinement interval
	    find_indx = find(abscissa >= x0);
	    indx = min(find_indx);
% 	    tmp = [[NaN*zeros(last_n,1);abscissa(1:indx-1)], ...
% 			[NaN*zeros(last_n,1);ordinate(1:indx-1)]];
% 	    last_few = tmp(indx:indx+last_n-1,:);
	    indx = indx - 1;   % it gets incremented below
	    % if calculating points is expensive, don't throw them away
	    if options.expensive_function
	      saved_points(saved_n+1:saved_n+length(find_indx),:) = ...
		  [abscissa(find_indx),ordinate(find_indx)];
	      saved_n = saved_n + length(find_indx);
	    end
	  end
	end
      end
    end

    % if at end of a refinement interval, drop down a level
    refine_indx(refine_level) = refine_indx(refine_level) + 1;
    while refine_indx(refine_level) > refine_steps(refine_level)
      refine_level = refine_level - 1;
      if refine_level == 0 break; end
      refine_indx(refine_level) = refine_indx(refine_level) + 1;
    end

    indx = indx + 1;

  end


% if an error occurs, output all successfully calculated points to
% save losing them, and set error flag if function call requested
% it
catch
  exitflag = -indx;
  % indx might be screwed so ensure last point is dumped to output
  indx = indx + 1;
  disp('Fatal error in SMARTPLOT:');
  disp(lasterr);
  disp('Dumping calculated points and aborting..');
end


% sort points by ascending abscissa value (otherwise
% plot(abscissa,ordinate) might be a mess!)
abscissa = abscissa(1:indx-1);
ordinate = ordinate(1:indx-1);
[abscissa,indx] = sort(abscissa);
ordinate = ordinate(indx);

varargout{1} = abscissa;
varargout{2} = ordinate;
if nargout == 3
  varargout{3} = exitflag;
end

