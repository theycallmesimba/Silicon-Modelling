classdef twoDimLOHO < handle
    %TWODIMWF Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        wavefunctionMG; % Wavefunction amplitude (in mesh grid format)
        wavefunctionNO; % Wavefunction amplitude (in natural ordering format)
        n; % nth mode in x axis
        m; % mth mode in y axis
        energy; % Wavefunction energy (en + em)
        dot; % Index corresponding to which dot wf is centered on
    end
    
    methods
        function obj = initialize(obj,wMG,wNO,nn,mm,en,di)
            obj.wavefunctionMG = wMG;
            obj.wavefunctionNO = wNO;
            obj.n = nn;
            obj.m = mm;
            obj.energy = en;
            obj.dot = di;
        end
    end
    
end

