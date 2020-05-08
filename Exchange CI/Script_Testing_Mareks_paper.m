A = [0.5750, 0.5750, 0.5750; 0.7086, -0.7086, 0; 0.4091, 0.4091, -0.8182];

nItinOrb = 3;
for jj = 1:nItinOrb
    for ii = 1:nItinOrb
        for kk = 1:nItinOrb
            for ll = 1:nItinOrb
                fprintf(1,'<%d%d|%d%d>(%d,%d)  ',ii,jj,kk,ll,(jj-1)*nItinOrb + ii, (kk-1)*nItinOrb + ll);
            end
        end
        fprintf(1,'\n');
    end
end
fprintf(1,'\n');
%%
t1 = 2.939179;
t2 = 0.512857;
t3 = 0.004653;
t4 = 0.003446;
t5 = 0.000049;
t6 = 0.000019;

CLOHO = zeros(9);
% On-site repulsion <rr|v|rr>
CLOHO(1,1) = t1; CLOHO(5,5) = t1; CLOHO(9,9) = t1;
% Inter-dot repulsion <rs|sr>
CLOHO(2,2) = t2; CLOHO(3,3) = t2; CLOHO(4,4) = t2; 
CLOHO(6,6) = t2; CLOHO(7,7) = t2; CLOHO(8,8) = t2;
% Small electron scattering terms
% Type 1 <rr|v|rs>
CLOHO(1,2) = t3; CLOHO(1,3) = t3;
CLOHO(5,4) = t3; CLOHO(5,6) = t3;
CLOHO(9,7) = t3; CLOHO(9,8) = t3;
% Type 2 <rs|v|st>
CLOHO(2,3) = t4; CLOHO(4,6) = t4; CLOHO(7,8) = t4;
% Type 3 <rs|v|rs>
CLOHO(2,4) = t5; CLOHO(3,7) = t5; CLOHO(6,8) = t5;
% Type 4 <rr|v|st>
CLOHO(1,6) = t6; CLOHO(1,8) = t6;
CLOHO(5,3) = t6; CLOHO(5,7) = t6;
CLOHO(9,2) = t6; CLOHO(9,3) = t6;
% Type 1* <rr|v|sr>
% CLOHO(1,4) = t3; CLOHO(1,7) = t3;  
% CLOHO(5,2) = t3; CLOHO(5,8) = t3;
% CLOHO(9,3) = t3; CLOHO(9,6) = t3;
% Type 4* <rr|v|ss>
% CLOHO(1,5) = t6; CLOHO(1,9) = t6; CLOHO(5,9) = t6;

tempMatrix = CLOHO - diag(diag(CLOHO));
CLOHO = CLOHO + tempMatrix';
for ii = 1:3
    for jj = 1:3
        rowInd = jj + (ii-1)*3;
        for kk = 1:3
            for ll = 1:3
                colInd = ll + (kk-1)*3;
                fprintf(1,'%.6f ',CLOHO(rowInd, colInd));
            end
        end
        fprintf(1,'\n');
    end
end

sparams.numElectrons = 3;
sparams.nSingleOrbitals = 3;
sparams.nItinerantOrbitals = sparams.nSingleOrbitals;
sparams.nOriginHOs = sparams.nItinerantOrbitals;

sparams.LCoriginEns = [-8.713,-8.6179,-8.6179];
% sparams.CMEsItin = kron(A,A)*CLOHO*kron(inv(A),inv(A));
% sparams.CMEsItin = kron(A,A)*CLOHO*kron(A',A');
% sparams.CMEsItin = kron(inv(A),inv(A))*CLOHO*kron(A,A);
%%
sparams.numElectrons = 3;
sparams.nSingleOrbitals = 3;
sparams.nItinerantOrbitals = sparams.nSingleOrbitals;
sparams.nOriginHOs = sparams.nItinerantOrbitals;

sparams.LCoriginEns = [-8.713,-8.6179,-8.6179];
sparams.CMEsItin = kron(A,A)*squeeze(storedCMEsLOHO(end,:,:))*kron(A',A');
size(sparams.CMEsItin)
% sparams.CMEsItin = squeeze(storedCMEsItin(end,:,:));
sparams = buildSecondQuantizationHam(sparams,[2]);

[eVectors, ens] = eigs(sparams.H2ndQ,8,'sa');
eVectors
ens
