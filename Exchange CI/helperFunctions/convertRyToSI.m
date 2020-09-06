function [conv_energy, conv_dist] = convertRyToSI(sparams, energy, dist)
%CONVERTUNITS Summary of this function goes here
%   Detailed explanation goes here
    hbar = 6.626E-34/2/pi;
    if strcmp(sparams.materialSystem,'Silicon')
        me = 9.10938356E-31*0.191;
    elseif strcmp(sparams.materialSystem,'GaAs')
        me = 9.10938356E-31*0.067; 
    end

    effaB = 4*pi*sparams.eps*(hbar^2)/(me*(sparams.ee^2)); % [m]
    effRy = (sparams.ee^2)/(8*pi*sparams.eps*effaB); % [J]

    conv_energy = energy*effRy;
    conv_dist = dist*effaB;
end

