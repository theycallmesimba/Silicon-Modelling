function effHamiltonianParams = buildEffHamiltonianParamVariable(epsL, epsR,...
    tc, Ez, Ex, deltaL, deltaR, S1, S2, EL, ER)
%BUILDEFFHAMILTONIANPARAMVARIABLE Summary of this function goes here
%   Detailed explanation goes here

    % TODO: There is a weird bug when the spin orbit is 0 that causes the
    % adiabatic parameter to have odd jumps and discontinuities.  We avoid
    % this by making the spin orbit strength, very small, but not exactly
    % 0.
    if S1 == 0
        S1 = 1.602E-31; % [10 aeV];
    end
    if S2 == 0
        S2 = 1.602E-31; % [10 aeV];
    end

    effHamiltonianParams{1} = epsL;
    effHamiltonianParams{2} = epsR;
    effHamiltonianParams{3} = tc;
    effHamiltonianParams{4} = Ez;
    effHamiltonianParams{5} = Ex;  
    effHamiltonianParams{6} = deltaL;
    effHamiltonianParams{7} = deltaR;
    effHamiltonianParams{8} = S1;
    effHamiltonianParams{9} = S2;
    effHamiltonianParams{10} = EL;
    effHamiltonianParams{11} = ER;
end

