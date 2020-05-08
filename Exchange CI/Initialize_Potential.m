fprintf(1,'************************\n');
fprintf(1,'* Loading potential... *\n');
fprintf(1,'************************\n\n');

sparams.unitsType = 'Rydberg';
% sparams.unitsType = 'SI - J';
% sparams.unitsType = 'SI - meV';
simparams;

if strcmp(sparams.unitsType,'Rydberg')
    effaB = 1; % [aB*]
    effRy = 1; % [Ry*]
% elseif strcmp(sparams.unitsType,'SI - meV')
%     effaB = 4*pi*sparams.eps*(sparams.hbar^2)/(sparams.me*(sparams.ee^2)); % [m]
%     effRy = (sparams.ee^2)/(8*pi*sparams.eps*effaB); % [meV]
elseif strcmp(sparams.unitsType,'SI - J')
    effaB = 4*pi*sparams.eps*(sparams.hbar^2)/(sparams.me*(sparams.ee^2)); % [m]
    effRy = (sparams.ee^2)/(8*pi*sparams.eps*effaB); % [J]
end

% Scale the dot locations according to our units
sparams.dotLocations = sparams.dotLocations*effaB;

% Make the potential wells
Vimin = 10;
Vimin = Vimin*effRy;
di = 2.3*effaB;
if strcmp(sparams.unitsType,'Rydberg')
    omega = sqrt(Vimin/di^2);
elseif strcmp(sparams.unitsType,'SI - meV') || strcmp(sparams.unitsType,'SI - J')
    omega = sqrt(2*Vimin/(sparams.me*di^2));
end

% Fill in grid parameters
gparams.ngridx = 250;
gparams.ngridy = 250;
gparams.xx = linspace(-10,10,gparams.ngridx)*effaB;
gparams.yy = linspace(-3,3,gparams.ngridy)*effaB;
[gparams.XX,gparams.YY] = meshgrid(gparams.xx,gparams.yy);

gparams.VV = zeros(gparams.ngridy,gparams.ngridx);
for ii = 1:sparams.nDots
    gparams.VV = gparams.VV + -Vimin*exp(-((gparams.XX - sparams.dotLocations(ii,1)).^2 +...
        (gparams.YY - sparams.dotLocations(ii,2)).^2)/di^2);
    sparams.fittedPotentialParameters(ii,:) = [omega, Vimin, sparams.dotLocations(ii,:)];
end

debugHere = 0;
if debugHere == 1
    plotMeshgrid(gparams,gparams.VV);
    set(gca,'visible','off')
end
fprintf(1,'Done!\n');