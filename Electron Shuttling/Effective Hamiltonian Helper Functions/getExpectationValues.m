function [orbExp, valExp, spinExp] = getExpectationValues(sparams,rho)
%GETEXPECTATIONVALUES Summary of this function goes here
%   Detailed explanation goes here
    orbExp = zeros(3,1);
    valExp = zeros(3,1);
    spinExp = zeros(3,1);
    
    if ~sparams.includeExcitedOrbital
        if sparams.includeOrbital
            if sparams.includeValley
                if sparams.includeSpin
                    orbExp(1,1) = trace(rho*kron(sparams.sigmax,eye(4)));
                    orbExp(2,1) = trace(rho*kron(sparams.sigmay,eye(4)));
                    orbExp(3,1) = trace(rho*kron(sparams.sigmaz,eye(4)));
                    valExp(1,1) = trace(rho*kron(eye(2),kron(sparams.sigmax,eye(2))));
                    valExp(2,1) = trace(rho*kron(eye(2),kron(sparams.sigmay,eye(2))));
                    valExp(3,1) = trace(rho*kron(eye(2),kron(sparams.sigmaz,eye(2))));
                    spinExp(1,1) = trace(rho*kron(eye(4),sparams.sigmax));
                    spinExp(2,1) = trace(rho*kron(eye(4),sparams.sigmay));
                    spinExp(3,1) = trace(rho*kron(eye(4),sparams.sigmaz));
                else
                    orbExp(1,1) = trace(rho*kron(sparams.sigmax,eye(2)));
                    orbExp(2,1) = trace(rho*kron(sparams.sigmay,eye(2)));
                    orbExp(3,1) = trace(rho*kron(sparams.sigmaz,eye(2)));
                    valExp(1,1) = trace(rho*kron(eye(2),sparams.sigmax));
                    valExp(2,1) = trace(rho*kron(eye(2),sparams.sigmax));
                    valExp(3,1) = trace(rho*kron(eye(2),sparams.sigmax));
                end
            elseif sparams.includeSpin
                orbExp(1,1) = trace(rho*kron(sparams.sigmax,eye(2)));
                orbExp(2,1) = trace(rho*kron(sparams.sigmay,eye(2)));
                orbExp(3,1) = trace(rho*kron(sparams.sigmaz,eye(2)));
                spinExp(1,1) = trace(rho*kron(eye(2),sparams.sigmax));
                spinExp(2,1) = trace(rho*kron(eye(2),sparams.sigmax));
                spinExp(3,1) = trace(rho*kron(eye(2),sparams.sigmax));
            else
                orbExp(1,1) = trace(rho*sparams.sigmax);
                orbExp(2,1) = trace(rho*sparams.sigmay);
                orbExp(3,1) = trace(rho*sparams.sigmaz);
            end
        elseif sparams.includeValley
            if sparams.includeSpin
                valExp(1,1) = trace(rho*kron(sparams.sigmax,eye(2)));
                valExp(2,1) = trace(rho*kron(sparams.sigmay,eye(2)));
                valExp(3,1) = trace(rho*kron(sparams.sigmaz,eye(2)));
                spinExp(1,1) = trace(rho*kron(eye(2),sparams.sigmax));
                spinExp(2,1) = trace(rho*kron(eye(2),sparams.sigmax));
                spinExp(3,1) = trace(rho*kron(eye(2),sparams.sigmax));
            else
                valExp(1,1) = trace(rho*sparams.sigmax);
                valExp(2,1) = trace(rho*sparams.sigmay);
                valExp(3,1) = trace(rho*sparams.sigmaz);
            end
        elseif sparams.includeSpin
            spinExp(1,1) = trace(rho*sparams.sigmax);
            spinExp(2,1) = trace(rho*sparams.sigmay);
            spinExp(3,1) = trace(rho*sparams.sigmaz);
        end
    end
end

