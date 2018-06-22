function [ xx, pots ] = loadPotentialsFromFigure( sparams )
%LOADPOTENTIALSFROMFIGURE Summary of this function goes here
%   Detailed explanation goes here

    nn = 1;
    for ii = 1:4
        for jj = 0:50
            currPotString = ['P' num2str(ii) '_' num2str(jj) '.fig'];
             
            % Load the figure
            fig = open([sparams.potDir currPotString]);
            axes = get(fig,'Children');
            data = get(axes,'Children');
            xdata = get(data,'XData');
            ydata = get(data,'YData');
            zdata = get(data,'ZData');
                        
            % Then we are on the first potential
            if ii == 1 && jj == 0
                unsampledxx = xdata;
                unsampledxx = unsampledxx';
                zz = ydata;
                
                % Get z position closest to 0
                [~,index] = min(abs(zz - 0));
                
                % The data may not be uniform in sampling, so we need to fix that for
                % the fourier transforms in the main code for speed up.
                % Find next highest power of 2 to the length of xx
                desiredGrid = 2^(nextpow2(length(unsampledxx)));
                % Make linearly spaced grid of points to interpolate the
                % potential at
                xx = linspace(min(xdata),max(xdata),desiredGrid);
                
                pots = zeros(4*51,length(xx));
            end
            pots(nn,:) = interp1(xdata,-zdata(index,:)*sparams.ee,xx);
            nn = nn + 1;
            close(fig);
        end
    end
    xx = xx*1E-9; % Convert to m
end

