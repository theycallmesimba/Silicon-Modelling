function sparams = solveCMEsSameOrbital( sparams )
%SOLVECMESSAMEORBITAL Summary of this function goes here
%   Detailed explanation goes here

    % First step is to find all the theta integrals first as they will be
    % non-zero for even power of aa  
    % Worst case for a on cos term is that na,nb,ng,nd = nNonShiftHOs with p's = 0
    % Likewise for b and m's.
    thetaMaxIndex = 4*(sparams.nNonShiftedHOs+1); % +1 for 0 cases 
    
    h = waitbar(0,sprintf('Theta indices a:%04d/%d  b:%04d/%d',0,thetaMaxIndex,0,thetaMaxIndex),...
        'Name','Finding theta integrals...',...
        'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
    
    sparams.thetaIntegrals = zeros(thetaMaxIndex+1,thetaMaxIndex+1);

    flag = 0;
    for aa = linspace(0,thetaMaxIndex,thetaMaxIndex+1)
        % Check if a is even, otherwise the integral is 0 
        if mod(aa,2) == 0
            for bb = linspace(0,thetaMaxIndex,thetaMaxIndex+1)
                %Check for cancel button click
                if getappdata(h,'canceling')
                    flag = 1;
                    break;
                end

                % Check if b is even, otherwise the integral is 0 
                if mod(bb,2) == 0
                    fun = @(th,a,b) (cos(th).^a).*(sin(th).^b);
                    temp = integral(@(theta)fun(theta,aa,bb),0,2*pi);
                    sparams.thetaIntegrals(aa+1,bb+1) = temp;
                end
                % If the integral is less than threshold, then we will consider
                % that to be 0 and can skip the rest of the integrals as
                % they will all be lower in magnitude.
                if temp < sparams.thetaIntThreshold
                    break;
                end
            end
            % Update waitbar
            waitbar(aa/thetaMaxIndex, h,...
                sprintf('Theta indices a:%04d/%d  b:%04d/%d',aa,thetaMaxIndex,bb,thetaMaxIndex));
        end
        if flag == 1
            break;
        end;
    end
    delete(h);
    
    % We need to loop through all of the possible two electron
    % configurations.  <alpha,beta|v|gamma,delta>
    % This part will be very long...
    sparams.CMEsNonShifted = zeros(sparams.nNonShiftedHOs^2,sparams.nNonShiftedHOs^2);
    
    % Produce lookup table for n and m values depending on what HO index we
    % are on.  This is needed for parallel computation in which we'll get
    % rid of any reference to sparams here since we don't want the huge
    % overhead in sending sparams to each worker in our parallel pool.
    indexLookup = zeros(sparams.nNonShiftedHOs,2);
    for ii = 1:sparams.nNonShiftedHOs
        indexLookup(ii,1) = sparams.nonShiftedHOs(ii).n;
        indexLookup(ii,2) = sparams.nonShiftedHOs(ii).m;
    end
    thetaInts = sparams.thetaIntegrals;
    numOfWFs = sparams.nNonShiftedHOs;
    
    % Calculate the scaling parameter
    omega = abs(mean(sparams.fittedPotentialParameters(:,1)));
    A = sqrt(sparams.hbar/(sparams.me*omega));
    
    h = waitbar(0,sprintf('alpha:%03d/%d  delta:%03d/%d',0,numOfWFs,0,numOfWFs),...
        'Name','Finding non shifted HO CMEs...',...
        'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
    ii = 0;
    
    for alpha = 1:numOfWFs
        nalpha = indexLookup(alpha,1);
        malpha = indexLookup(alpha,2);
        for beta = 1:numOfWFs
            ii = ii + 1;
            
            % Update waitbar
            waitbar(ii/(numOfWFs^2), h,...
                sprintf('alpha:%03d/%d  beta:%03d/%d',alpha,numOfWFs,beta,numOfWFs));
            
            nbeta = indexLookup(beta,1);
            mbeta = indexLookup(beta,2);
                        
            %Check for cancel button click
            if getappdata(h,'canceling')
                flag = 1;
                break;
            end
            
            tmpGammaDelta = zeros(numOfWFs,numOfWFs);
            
            parfor gamma = 1:numOfWFs  
%             for gamma = 1:numOfWFs  
                ngamma = indexLookup(gamma,1);
                mgamma = indexLookup(gamma,2);
                
                tmpDelta = zeros(1,numOfWFs);
                
                for delta = 1:numOfWFs
                    ndelta = indexLookup(delta,1);
                    mdelta = indexLookup(delta,2);
                    
%                     if alpha == 3 && beta == 3 && gamma == 3 && delta == 2
%                     fprintf(1,'\n%d,%d  %d,%d  %d,%d  %d,%d\n',...
%                         nalpha,malpha,nbeta,mbeta,ngamma,mgamma,ndelta,mdelta);
                    
                    % Do the embedded loops
                    result = CMEhelper(thetaInts,nalpha,malpha,nbeta,mbeta,...
                        ngamma,mgamma,ndelta,mdelta);
                    
                    result = result/(2*pi*A*sqrt(2));
                    result = result/sqrt(factorial(nalpha)*factorial(malpha)*...
                        factorial(nbeta)*factorial(mbeta)*factorial(ngamma)*...
                        factorial(mgamma)*factorial(ndelta)*factorial(mdelta));
%                     fprintf(1,'CME value = %f\n',result);
%                         close all force;
%                         return;
%                     end
%                     result = 0;
                    % Making a row vector with elements increasing as
                    % [<a_0b_0|v|g_0d_1> ... <a_0b_0|v|g_0d_M>]
                    % M = sparams.nNonShiftedHOs
                    tmpDelta(delta) = result;    
                end
                % Matrix with elements
                % [<a_0b_0|v|g_1d_1>, <a_0b_0|v|g_1d_2> ... <a_0b_0|v|g_1d_M>]
                % [<a_0b_0|v|g_2d_1>, <a_0b_0|v|g_2d_2> ... <a_0b_0|v|g_2d_M>]
                % [       ...                ...                   ...       ]
                % [<a_0b_0|v|g_Md_1>, <a_0b_0|v|g_Md_2> ... <a_0b_0|v|g_Md_M>]
                % M = sparams.nNonShiftedHOs
                tmpGammaDelta(gamma,:) = tmpDelta;
            end

            % Convert to a row vector with elements increasing as
            % [<a_0b_0|v|g_1d_1> ... <a_0b_0|v|g_1d_M>, <a_0b_0|v|g_2d_1> ... <a_0b_0|v|g_Md_M>]
            % M = sparams.nNonShiftedHOs
            temp = reshape(tmpGammaDelta',[1,numOfWFs^2]);
            
            % Matrix with elements
            % [<a_1b_1|v|g_1d_1>, <a_1b_1|v|g_1d_2> ... <a_1b_1|v|g_1d_M>, <a_1b_1|v|g_2d_1> ... <a_1b_1|v|g_Md_M>]
            % [<a_1b_2|v|g_1d_1>, <a_1b_2|v|g_1d_2> ... <a_1b_2|v|g_1d_M>, <a_1b_2|v|g_2d_1> ... <a_1b_2|v|g_Md_M>]
            % [       ...                ...                   ...                ...                   ...       ]
            % [<a_1b_M|v|g_1d_1>, <a_1b_M|v|g_1d_2> ... <a_1b_M|v|g_1d_M>, <a_1b_M|v|g_2d_1> ... <a_1b_M|v|g_Md_M>]
            % [<a_2b_1|v|g_1d_1>, <a_2b_1|v|g_1d_2> ... <a_2b_1|v|g_1d_M>, <a_2b_1|v|g_2d_1> ... <a_2b_1|v|g_Md_M>]
            % [       ...                ...                   ...                ...                   ...       ]
            % [<a_Mb_M|v|g_1d_1>, <a_Mb_M|v|g_1d_2> ... <a_Mb_M|v|g_1d_M>, <a_Mb_M|v|g_2d_1> ... <a_Mb_M|v|g_Md_M>]
            % M = sparams.nNonShiftedHOs
            sparams.CMEsNonShifted(beta + (alpha - 1)*numOfWFs,:) = temp;
            
            if flag == 1
                break;
            end;
        end
        if flag == 1
            break;
        end;
    end
    delete(h);
    
    % As a final step, we need to multiply our CME by the scaling factor from
    % our Coulomb potential e^2/(4*pi*e_r)
    k = 1/(4*pi*sparams.eps)*sparams.ee*sparams.ee;
    sparams.CMEsNonShifted = k*sparams.CMEsNonShifted;
end

% Just a function to help reduce clutter up above with all the summations
% going on
function temp = CMEhelper(thetaIntegrals, na, ma, nb, mb, ng, mg, nd, md)

    temp = 0;
    
    % Here we do a nice check to skip a bunch of potential loops.  Since
    % the p's are multiplied by 2's in the aa and bb calculation, the only
    % way the theta integrals will be non-zero are if na + nd + nb + ng is
    % even AND ma + md + mb + mg is even.  If these conditions aren't met,
    % then all of the theta integrals will be 0, so the CME element is 0.
    if mod(na + nd + nb + ng,2) ~= 0 || mod(ma + md + mb + mg,2) ~= 0
        temp = 0;
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
                    
%                     fprintf('aa = %d, bb = %d\n',aa,bb);
                    
                    % Look up theta integral table and see if we can skip
                    % this or not
                    tempTheta = thetaIntegrals(aa+1,bb+1);
%                     if tempTheta == 0
%                         continue;
%                     end
                    coef4 = coef3*factorial(p4)*nchoosek(mb,p4)*nchoosek(mg,p4);

                    p = (aa + bb)/2;
                    temp = temp + coef4*tempTheta*gamma(p + 1/2);      
                end
            end
        end
    end
end



