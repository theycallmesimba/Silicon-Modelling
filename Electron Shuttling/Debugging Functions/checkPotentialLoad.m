function checkPotentialLoad( sparams, xx, zz )

%CHECKPOTENTIALLOAD Summary of this function goes here
%   Detailed explanation goes here
    test = 0.6:0.01:0.8;

    fig = figure;
    for ii = 1:length(test)
        clf;
        testPot = sparams.P2DEGInterpolant({0.6,test(ii),0.6,0.7,0.7,xx});
        testPot = squeezeFast(sparams.numOfGates,testPot);
        plot(xx,testPot/sparams.ee);
        pause(0.2);
    end
    delete(fig);
    
    [XX,ZZ] = meshgrid(xx,zz);
    fig = figure;
    testPot = sparams.P2DInterpolant({0.6,0.6,0.768,0.6,0.6,zz,xx});
    testPot = squeezeFast(sparams.numOfGates,testPot);
    s = surf(XX,ZZ,testPot);
    set(s,'edgecolor','none');
    view(2);
    pause(5);
    delete(fig);
end

