for ii = 1:sparams.nSingleOrbitals
    fprintf('A norm of row %d: %f\n', ii, norm(sparams.acoeffs(ii,:)));
end
for ii = 1:sparams.nSingleOrbitals
    fprintf('B norm of row %d: %f\n', ii, norm(sparams.bcoeffs(ii,:)));
end
C = sparams.acoeffs*sparams.bcoeffs;
for ii = 1:sparams.nSingleOrbitals
    fprintf('C norm of row %d: %f\n', ii, norm(C(ii,:)));
end
Bid = sparams.bcoeffs*sparams.bcoeffs';
Aid = sparams.acoeffs*sparams.acoeffs';