classdef twoDimLCHO < handle
    %TWODIMWF Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        wavefunctionMG; % Wavefunction amplitude (in mesh grid format)
        wavefunctionNO; % Wavefunction amplitude (in mesh grid format)
        energy; % Wavefunction energy
    end
    
    methods
        function obj = initialize(obj,wMG,wNO,en)
            obj.wavefunctionMG = wMG;
            obj.wavefunctionNO = wNO;
            obj.energy = en;
        end
    end
    
end
