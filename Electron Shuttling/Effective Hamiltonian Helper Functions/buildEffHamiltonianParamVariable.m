function effHamiltonianParams = buildEffHamiltonianParamVariable(epsL, epsR,...
    tc, Ez, Ex, deltaL, deltaR, S1, S2)
%BUILDEFFHAMILTONIANPARAMVARIABLE Summary of this function goes here
%   Detailed explanation goes here

    effHamiltonianParams{1} = epsL;
    effHamiltonianParams{2} = epsR;
    effHamiltonianParams{3} = tc;
    effHamiltonianParams{4} = Ez;
    effHamiltonianParams{5} = Ex;  
    effHamiltonianParams{6} = deltaL;
    effHamiltonianParams{7} = deltaR;
    effHamiltonianParams{8} = S1;
    effHamiltonianParams{9} = S2;
end

