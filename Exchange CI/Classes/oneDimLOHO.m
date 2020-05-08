classdef oneDimLOHO < handle
    %LOCAL1DWF Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        wavefunction; % Wavefunction amplitude vector
        n; % nth mode
        energy; % Wavefunction energy
        dot; % Index corresponding to which dot wf is centered on
    end
    
    methods
        function obj = initialize(obj,w,nn,en,di)
            obj.wavefunction = w;
            obj.n = nn;
            obj.energy = en;
            obj.dot = di;
        end
    end
    
end

