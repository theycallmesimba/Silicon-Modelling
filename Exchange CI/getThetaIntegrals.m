function thetaIntegrals = getThetaIntegrals(sparams)
%GETTHETAINTEGRALS Summary of this function goes here
%   Detailed explanation goes here

    % First step is to find all the theta integrals first as they will be
    % non-zero for even power of aa  
    % Worst case for a on cos term is that na,nb,ng,nd = nNonShiftHOs with p's = 0
    % Likewise for b and m's.
    thetaMaxIndex = 4*(sparams.nOriginHOs+1); % +1 for 0 cases 
    
    h = waitbar(0,sprintf('Theta indices a:%04d/%d  b:%04d/%d',0,thetaMaxIndex,0,thetaMaxIndex),...
        'Name','Finding theta integrals...',...
        'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
    
    thetaIntegrals = zeros(thetaMaxIndex+1,thetaMaxIndex+1);

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
                    thetaIntegrals(aa+1,bb+1) = temp;
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
        end
    end
    delete(h);
end



