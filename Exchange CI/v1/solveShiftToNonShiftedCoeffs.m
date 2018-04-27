function [ sparams ] = solveShiftToNonShiftedCoeffs( sparams )
%SOLVESHIFTTONONSHIFTEDCOEFFS Summary of this function goes here
%   Detailed explanation goes here    
    sparams.bcoeffs = zeros(sparams.nSingleOrbitals,...
        (sparams.maxNonShiftedHOsX+1)*(sparams.maxNonShiftedHOsY+1));
    
    h = waitbar(0,'1','Name','Unfolding Shifted HOs...',...
        'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
    
    kk = 0;
    flag = 0;
    for ii = 1:((sparams.maxNonShiftedHOsX+1)*(sparams.maxNonShiftedHOsY+1))
        nonshiftedNO = sparams.nonShiftedHOs(ii).wavefunctionNO;
        for jj = 1:sparams.nSingleOrbitals
            %Check for cancel button click
            if getappdata(h,'canceling')
                flag = 1;
                break;
            end
            shiftedNO = sparams.localHOs(jj).wavefunctionNO;
            
            % B matrix with B_ij = <r_j|alpha_i>
            sparams.bcoeffs(jj,ii) = sum(conj(shiftedNO).*nonshiftedNO);
            
            kk = kk + 1;
            % Update waitbar
            waitbar(kk/(sparams.nSingleOrbitals*(sparams.maxNonShiftedHOsX+1)*(sparams.maxNonShiftedHOsY+1)), h,...
                sprintf('Shifted Index:%04d/%d  Non-shifted Index:%04d/%d',jj,(sparams.maxNonShiftedHOsX+1)*(sparams.maxNonShiftedHOsY+1)),ii,sparams.nSingleOrbitals);
        end
        if flag == 1
            break;
        end
    end
    
    % Normalize each row now
    for ii = 1:sparams.nSingleOrbitals
        sparams.bcoeffs(ii,:) = sparams.bcoeffs(ii,:)/norm(sparams.bcoeffs(ii,:));
    end
    
    % Close waitbar
    delete(h);
end

