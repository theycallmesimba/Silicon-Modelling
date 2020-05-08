function CME = calculateOriginCME(na, ma, nb, mb, ng, mg, nd, md)

    CME = 0;
 
    % Intialize a and p
    a0 = na + nb + ng + nd;
    b = ma + mb + mg + md;
    pInd0 = a0 + ma + mb + mg + md;
    
    % If a and b are both even, then the CME will be non-zero.  Otherwise,
    % it will be zero and we can just skip doing all of these loops.  Note
    % that we don't need to consider the -2*p_i terms in a and b because
    % those are even.
    if mod(a0,2) ~= 0 || mod(b,2) ~= 0
        return
    end
    
    for p1 = 0:min(na,nd) 
        coef1 = factorial(p1)*nchoosek(na,p1)*nchoosek(nd,p1);
        % Continue building a and p
        a1 = a0 - 2*p1;
        pInd1 = pInd0 - 2*p1;
        
        for p2 = 0:min(ma,md)
            coef2 = coef1*factorial(p2)*nchoosek(ma,p2)*nchoosek(md,p2);
            % Continue building p
            pInd2 = pInd1 - 2*p2;
            
            for p3 = 0:min(nb,ng)
                coef3 = coef2*factorial(p3)*nchoosek(nb,p3)*nchoosek(ng,p3);
                % Finish building a and continue building and p
                a = a1 - 2*p3;
                pInd3 = pInd2 - 2*p3;
        
                for p4 = 0:min(mb,mg)                    
                    coef4 = coef3*factorial(p4)*nchoosek(mb,p4)*nchoosek(mg,p4);
                    % Finish building p
                    p = (pInd3 - 2*p4)/2;
                    
                    % Skip this sum term if 2p is odd
                    if floor(p) ~= p
                        continue
                    end
                    
                    % Calculate the CME
                    CME = CME + (-1)^p*coef4*gamma(p + 1/2)*...
                        beta(p - ((a - 1)/2),(a + 1)/2); 
                end
            end
        end
    end
    
    % Take care of the front coefficients
    globPhase = (-1)^(nb + mb + ng + mg);
    CME = CME*2/(pi*sqrt(2))*globPhase;
    CME = CME/sqrt(factorial(na)*factorial(ma)*...
        factorial(nb)*factorial(mb)*factorial(ng)*...
        factorial(mg)*factorial(nd)*factorial(md));
end