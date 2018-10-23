function [epsL, epsR, tc, Ez, Ex, deltaL, deltaR, S1, S2, EL, ER] =...
    decodeEffHamiltonianParamVariable(effHamiltonianParams)
%DECODEEFFHAMILTONIANPARAMVARIABLE Summary of this function goes here
%   Detailed explanation goes here

    epsL = effHamiltonianParams{1};
    epsR = effHamiltonianParams{2};
    tc = effHamiltonianParams{3};
    Ez = effHamiltonianParams{4};
    Ex = effHamiltonianParams{5}; 
    deltaL = effHamiltonianParams{6};
    deltaR = effHamiltonianParams{7};
    S1 = effHamiltonianParams{8};
    S2 = effHamiltonianParams{9};
    EL = effHamiltonianParams{10};
    ER = effHamiltonianParams{11};
end

