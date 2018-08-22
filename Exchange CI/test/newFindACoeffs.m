function sparams = newFindACoeffs(sparams,X,Y,V)
%NEWFINDACOEFFS Summary of this function goes here
%   Detailed explanation goes here

    [Rwfs, ~] = solve2DSingleElectronSE(sparams, X, Y, V, sparams.nNonShiftedHOs);
    
    sparams.acoeffs1 = zeros(sparams.nNonShiftedHOs);
    for ii = 1:sparams.nNonShiftedHOs
        for jj = 1:sparams.nNonShiftedHOs
%             size(sparams.nonShiftedHOs(ii).wavefunctionMG)
%             size(Rwfs(jj,:,:))
            sparams.acoeffs1(ii,jj) = getInnerProduct(sparams.nonShiftedHOs(ii).wavefunctionMG,...
                Rwfs(:,:,jj),X,Y);
        end
    end
    for ii = 1:sparams.nNonShiftedHOs
        sparams.acoeffs1(ii,:) = sparams.acoeffs1(ii,:)./norm(sparams.acoeffs1(ii,:));
    end
    sparams.acoeffs1 = sparams.acoeffs1';;
end

