function Heff = constructEffectiveHamiltonian( sparams, epsL, epsR,...
    tc, Ez, deltaL, deltaR, S1, S2)
%CONSTRUCTEFFECTIVEHAMILTONIAN Summary of this function goes here
%   Detailed explanation goes here

    if sparams.includeOrbital
        Hdet = [epsL 0;0 epsR];
        Htun = [0 tc;tc 0];
        Horbital = Hdet + Htun;
        
        if sparams.includeValley
            Horbital = kron(Horbital,eye(2));
            
            Hvalley = kron([1 0;0 0],[0 deltaL;deltaL' 0]) + kron([0 0;0 1],[0 deltaR;deltaR' 0]);
            
            if sparams.includeSpin
                Horbital = kron(Horbital,eye(2));
                Hvalley = kron(Hvalley,eye(2));
                
                Hzeeman = [Ez 0;0 -Ez];
                Hzeeman = kron(eye(4),Hzeeman);
                
                projL = kron([1 0;0 0],eye(2)); projR = kron([0 0;0 1],eye(2));
                HspinOrbit1 = kron(projL,[0 S1;S1 0]) + kron(projR,[0 -S1;-S1 0]);
                HspinOrbit2 = kron(kron([0 1;0 0],eye(2)),[0 S2;0 0]) +...
                    kron(kron([0 0;1 0],eye(2)),[0 0;S2 0]) +...
                    kron(kron([0 0;1 0],eye(2)),[0 -S2;0 0]) +...
                    kron(kron([0 1;0 0],eye(2)),[0 0;-S2 0]);
                
                Heff = Horbital + Hvalley + Hzeeman + HspinOrbit1 + HspinOrbit2;
            end
        elseif sparams.includeSpin
            Horbital = kron(Horbital,eye(2));

            Hzeeman = [Ez 0;0 -Ez];
            Hzeeman = kron(eye(2),Hzeeman);
            
            projL = [1 0;0 0]; projR = [0 0;0 1];
            HspinOrbit1 = kron(projL,[0 S1;S1 0]) + kron(projR,[0 -S1;-S1 0]);
            HspinOrbit2 = kron([0 1;0 0],[0 S2;0 0]) + kron([0 0;1 0],[0 0;S2 0]) +...
                kron([0 0;1 0],[0 -S2;0 0]) + kron([0 1;0 0],[0 0;-S2 0]);
            
            Heff = Horbital + Hzeeman + HspinOrbit1 + HspinOrbit2;
        end
    elseif sparams.includeValley
        Hvalley = [deltaL 0;0 -deltaL];
        
        if sparams.includeSpin
            Hvalley = kron(Hvalley,eye(2));
            
            Hzeeman = [Ez 0;0 -Ez];
            Hzeeman = kron(eye(2),Hzeeman);
            
            Heff = Hvalley + Hzeeman;
        end
    elseif sparams.includeSpin
        Hzeeman = [Ez 0;0 -Ez];
        
        Heff = Hzeeman;
    end
end

