function [orbPurity, valPurity, spinPurity, totPurity] = getPurityValues(sparams, rho)
%GETPURITYVALUES Summary of this function goes here
%   Detailed explanation goes here
    orbPurity = 0;
    valPurity = 0;
    spinPurity = 0;
    totPurity = trace(rho^2);

    if ~sparams.includeExcitedOrbital
        if sparams.includeOrbital
            if sparams.includeValley
                if sparams.includeSpin
                    if sparams.includeSecondSpin
                        orbRho = partialTrace(rho,[2,3,4],[2,2,2,2]);
                        valRho = partialTrace(rho,[1,3,4],[2,2,2,2]);
                        spinRho = partialTrace(rho,[1,2],[2,2,2,2]);
                        
                        orbPurity = trace(orbRho^2);
                        valPurity = trace(valRho^2);
                        spinPurity = trace(spinRho^2);
                    else
                        orbRho = partialTrace(rho,[2,3],[2,2,2]);
                        valRho = partialTrace(rho,[1,3],[2,2,2]);
                        spinRho = partialTrace(rho,[1,2],[2,2,2]);
                        
                        orbPurity = trace(orbRho^2);
                        valPurity = trace(valRho^2);
                        spinPurity = trace(spinRho^2);
                    end
                else
                    orbRho = partialTrace(rho,[2],[2,2]);
                    valRho = partialTrace(rho,[1],[2,2]);

                    orbPurity = trace(orbRho^2);
                    valPurity = trace(valRho^2);
                end
            elseif sparams.includeSpin
                if sparams.includeSecondSpin
                    orbRho = partialTrace(rho,[2,3],[2,2,2]);
                    spinRho = partialTrace(rho,[1],[2,2,2]);

                    orbPurity = trace(orbRho^2);
                    spinPurity = trace(spinRho^2);
                else
                    orbRho = partialTrace(rho,[2],[2,2]);
                    spinRho = partialTrace(rho,[1],[2,2]);

                    orbPurity = trace(orbRho^2);
                    spinPurity = trace(spinRho^2);
                end
            else
                orbPurity = trace(rho^2);
            end
        elseif sparams.includeValley
            if sparams.includeSpin
                if sparams.includeSecondSpin
                    valRho = partialTrace(rho,[2,3],[2,2,2]);
                    spinRho = partialTrace(rho,[1],[2,2,2]);

                    valPurity = trace(valRho^2);
                    spinPurity = trace(spinRho^2);
                else
                    valRho = partialTrace(rho,[2],[2,2]);
                    spinRho = partialTrace(rho,[1],[2,2]);

                    valPurity = trace(valRho^2);
                    spinPurity = trace(spinRho^2);
                end
            else
                valPurity = trace(rho^2);
            end
        elseif sparams.includeSpin
            if sparams.includeSecondSpin
                spinPurity = trace(rho^2);               
            else
                spinPurity = trace(rho^2);
            end
        end
    end
end

