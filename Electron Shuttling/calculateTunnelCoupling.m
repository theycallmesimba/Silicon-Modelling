function tc = calculateTunnelCoupling( sparams, xx, vv )
%CALCULATETUNNELCOUPLING Summary of this function goes here
%   Detailed explanation goes here

    % We will use a semi-classical WKB approach to calculate the tunnel
    % coupling as laid out in https://arxiv.org/pdf/1205.0366.pdf
    
    % The first step is to find the tunnel barrier height and location
    [v0, locs] = findpeaks(vv);
%     findpeaks(vv);
%     locs
    if length(locs) > 1
%         fprintf(1,'ERROR: findpeaks found two peaks in the potential!\n');
        [v0, ind] = min(v0);
        locs = locs(ind);
%         tc = -1;
%         return;
    elseif length(locs) <= 0
        fprintf(1,'ERROR: findpeaks could not find any peaks in the potential!\n');
        tc = NaN;
        return;
    end

    % Now using the barrier as the "center" point.  Create two potentials
    % where one is the left side of the tunnel barrier mirrored about the
    % y-axis and one is the right side of the tunnel barrier mirrored about
    % the y-axis.
    potL = [vv(1:locs) fliplr(vv(1:locs-1))];
    [~,ens] = solve1DSingleElectronSE(sparams,2,[xx(1),xx(2)],potL);
    tcL = (ens(2,2) - ens(1,1))/2;
    potR = [fliplr(vv(locs+1:end)) vv(locs:end)];
    [~,ens] = solve1DSingleElectronSE(sparams,2,[xx(1),xx(2)],potR);
    tcR = (ens(2,2) - ens(1,1))/2;
    
    % Now find the minima in each of the wells
    [eps,~] = findpeaks(-vv);
    epsL = -eps(1);
    epsR = -eps(2);
    
    temp1 = ((v0 - epsL)/(v0 - epsR))^(1/4);
    temp2 = ((v0 - epsR)/(v0 - epsL))^(1/4);
    An = 1/2*(temp1 + temp2);
    
    tc = An*sqrt(tcL*tcR);
end

