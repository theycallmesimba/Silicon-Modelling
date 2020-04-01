function kets = setKetSameGlobalPhase(kets)
%SETKETSAMEGLOBALPHASE Summary of this function goes here
%   This funciton takes all of the kets in the argument and set them to all
%   have the same global phase where the first element of the kets is real
%   and positive.
    [~,nKets] = size(kets);
    
    for ii = 1:nKets
        if imag(kets(1,ii)) ~= 0 || real(kets(1,ii)) < 0 
            kets(:,ii) = kets(:,ii)*exp(-1i*angle(kets(1,ii)));
        end
    end
end

