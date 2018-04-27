classdef twoDimNonShiftHO < handle
    %TWODIMNONSHIFTHO Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        wavefunctionMG; % Wavefunction amplitude (in mesh grid format)
        wavefunctionNO; % Wavefunction amplitude (in natural ordering format)
        n; % nth mode in x axis
        m; % mth mode in y axis
    end
    
    methods
        function obj = initialize(obj,wMG,wNO,nn,mm)
            obj.wavefunctionMG = wMG;
            obj.wavefunctionNO = wNO;
            obj.n = nn;
            obj.m = mm;
        end
    end
    
end

