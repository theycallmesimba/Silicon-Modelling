function [ epsL, epsR ] = getDetuning( sparams, xx, vv, dotLocs )
%GETDETUNING Summary of this function goes here
    % We want to find the middle point that is defining the tunnel
    % coupling between the two dots.  This will let us create two
    % symmetric double well potentials (where we've flipped the potential
    % along the x-axis).  From this, we will find the tunnel coupling for
    % these potentials which should have 0 detuning.  Therefore, we can
    % extract what the local energy levels would be in the absence of
    % tunnel coupling

    % First though, we need to find which dots were are currently looking
    % at.  To do this, we'll look at the potential, find the two lowest
    % peak points and use those.
    [pks,locs] = findpeaks(-vv);
    [~,sortInd] = sort(-pks);
    locs = locs(sortInd);
        
    % First dot...
    [~,LDot] = min(abs(dotLocs - xx(locs(1))));
    % Second dot...
    [~,RDot] = min(abs(dotLocs - xx(locs(2))));
        
    midPoint = (dotLocs(LDot) + dotLocs(RDot))/2;
    [~,locOfInterest] = min(abs(xx - midPoint));
    
    dx = xx(2) - xx(1);
        
    % Get the symmetric left side potential
    potSymL = [vv(1:locOfInterest), fliplr(vv(1:(locOfInterest-1)))];
    xxSymL = ((1:length(potSymL))-1)*dx;
    % Get the symmetric right side potential
    potSymR = [fliplr(vv(locOfInterest:end)),vv(locOfInterest+1:end)];
    xxSymR = ((1:length(potSymR))-1)*dx;
    
    [~,ensSymL] = solve1DSingleElectronSE(sparams, 2, xxSymL, potSymL);
    [~,ensSymR] = solve1DSingleElectronSE(sparams, 2, xxSymR, potSymR);
    epsL = (ensSymL(1,1) + ensSymL(2,2))/2;
    epsR = (ensSymR(1,1) + ensSymR(2,2))/2;
end

