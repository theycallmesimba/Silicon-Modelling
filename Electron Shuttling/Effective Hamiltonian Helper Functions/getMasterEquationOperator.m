function dRho_dt = getMasterEquationOperator(sparams,rho,H,T2)
%GETMASTEREQUATIONOPERATOR Summary of this function goes here
%   Detailed explanation goes here
    unitaryTerm = -1i/sparams.hbar*(H*rho - rho*H);
    dRho_dt = unitaryTerm;
    
    linOp = sqrt(1/T2)*[0,0;0,1];
    if sparams.includeT2
        if sparams.includeValley
            linOp = kron(linOp,eye(2));
            if sparams.includeSpin
                linOp = kron(linOp,eye(2));
            end
        elseif sparams.includeSpin
            linOp = kron(linOp,eye(2));
        end
        
        dissipatorTerm = -1/2*((linOp')*linOp*rho + rho*(linOp')*linOp) + linOp*rho*linOp;
        dRho_dt = dRho_dt + dissipatorTerm;
    end
end


