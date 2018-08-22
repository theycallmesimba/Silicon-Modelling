function CME = calculateCME(thetaIntegrals, bohrRadii, k, na, ma, nb, mb, ng, mg, nd, md)

    CME = 0;
    
    % Here we do a nice check to skip a bunch of potential loops.  Since
    % the p's are multiplied by 2's in the aa and bb calculation, the only
    % way the theta integrals will be non-zero are if na + nd + nb + ng is
    % even AND ma + md + mb + mg is even.  If these conditions aren't met,
    % then all of the theta integrals will be 0, so the CME element is 0.
    if mod(na + nd + nb + ng,2) ~= 0 || mod(ma + md + mb + mg,2) ~= 0
        CME = 0;
        return;
    end
    
    for p1 = 0:min(na,nd) 
        coef1 = factorial(p1)*nchoosek(na,p1)*nchoosek(nd,p1);
        for p2 = 0:min(ma,md)
            coef2 = coef1*factorial(p2)*nchoosek(ma,p2)*nchoosek(md,p2);
            for p3 = 0:min(nb,ng)
                coef3 = coef2*factorial(p3)*nchoosek(nb,p3)*nchoosek(ng,p3);
                for p4 = 0:min(mb,mg)                     
                    aa = na + nd + nb + ng - 2*p1 - 2*p3;
                    bb = ma + md + mb + mg - 2*p2 - 2*p4;
                    
                    coef4 = coef3*factorial(p4)*nchoosek(mb,p4)*nchoosek(mg,p4);

                    p = (aa + bb)/2;
                    CME = CME + coef4*thetaIntegrals(aa+1,bb+1)*gamma(p + 1/2);      
                end
            end
        end
    end
    
    % Do some scalings of the terms
    CME = CME/(2*pi*bohrRadii*sqrt(2));
    CME = CME/sqrt(factorial(na)*factorial(ma)*...
        factorial(nb)*factorial(mb)*factorial(ng)*...
        factorial(mg)*factorial(nd)*factorial(md));
    CME = k*CME;
end