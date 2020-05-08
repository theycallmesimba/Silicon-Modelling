function updateShuttlingFigure( sparams, fig, wf1, wf2, vv )
%UPDATEFIGURE Summary of this function goes here
%   Detailed explanation goes here
    
    figure(fig);
    h = findobj('type','line');
    % Ground state WF
    h(1).YData = abs(wf2).^2/norm(abs(wf2).^2);
    % Simulated state WF
    h(2).YData = abs(wf1).^2/norm(abs(wf1).^2);
    % Potential
%     CBoffset = -2.2747;
%     vvTemp = vv/sparams.ee - CBoffset;
%     h(3).YData = vvTemp;
    h(3).YData = vv/sparams.ee;
end

