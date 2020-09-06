function CME = calculateOriginCMEOPTIMIZED(na, ma, nb, mb, ng, mg, nd, md,...
    factorialLookup, nchoosekLookup, gammaLookup, betaLookup)

    CME = 0;
 
    % Intialize a and 2p
    % Calculating a/2 instead of a to avoid repeated division in loops
    a0 = (na + nb + ng + nd);
    adiv2_0 = a0/2;
    b = ma + mb + mg + md;
    twoPInd0 = 2*adiv2_0 + b;
    % Calculating p instead of 2p to avoid repeated division in loops
    pInd0 = twoPInd0/2; 
    
    % If a and b are both even, then the CME will be non-zero.  Otherwise,
    % it will be zero and we can just skip doing all of these loops.  Note
    % that we don't need to consider the -2*p_i terms in a and b because
    % those are even.
    if mod(2*adiv2_0,2) ~= 0 || mod(b,2) ~= 0
        return
    end
    
    for p1 = 0:min(na,nd)
        coef1 = factorialLookup(p1+1)*nchoosekLookup(na+1,p1+1)*...
            nchoosekLookup(nd+1,p1+1);
        % Continue building a and p
        adiv2_1 = adiv2_0 - p1;
        pInd1 = pInd0 - p1;
        
        for p2 = 0:min(ma,md)
            coef2 = coef1*factorialLookup(p2+1)*nchoosekLookup(ma+1,p2+1)*...
                nchoosekLookup(md+1,p2+1);
            % Continue building p
            pInd2 = pInd1 - p2;
            
            for p3 = 0:min(nb,ng)
                coef3 = coef2*factorialLookup(p3+1)*nchoosekLookup(nb+1,p3+1)*...
                    nchoosekLookup(ng+1,p3+1);
                % Finish building a and continue building and p
                adiv2 = adiv2_1 - p3;
                pInd3 = pInd2 - p3;
        
                for p4 = 0:min(mb,mg)                    
                    coef4 = coef3*factorialLookup(p4+1)*nchoosekLookup(mb+1,p4+1)*...
                        nchoosekLookup(mg+1,p4+1);
                    % Finish building p
                    p = pInd3 - p4;
                    
                    % Calculate the CME
                    CME = CME + (-1)^p*coef4*gammaLookup(p+1)*...
                        betaLookup(p - adiv2 + 1, adiv2 + 1);
                end
            end
        end
    end
    
    % Take care of the front coefficients
    globPhase = (-1)^(nb + mb + ng + mg);
    CME = CME*2/(pi*sqrt(2))*globPhase;
    CME = CME/sqrt(factorialLookup(na+1)*factorialLookup(ma+1)*...
        factorialLookup(nb+1)*factorialLookup(mb+1)*factorialLookup(ng+1)*...
        factorialLookup(mg+1)*factorialLookup(nd+1)*factorialLookup(md+1));
end