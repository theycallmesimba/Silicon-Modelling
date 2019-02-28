function dPsi_dt = schrodingerEqODE(t, psi, sparams, xx)
    gPulse = getInterpolatedPulseValues(sparams,t,sparams.vPulseGInterpolants);
    vv = sparams.P2DEGInterpolant([num2cell(gPulse'),xx]);
    vv = squeezeFast(sparams.numOfGates,vv);
    currH = make1DSELap(sparams,xx,vv);
    
    dPsi_dt = -1i/sparams.hbar*currH*psi;
end