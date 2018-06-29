function interpArg = getInterpolantArgument( voltageVec, xx, zz )
%GETINTERPOLANTARGUMENT Summary of this function goes here
%   Detailed explanation goes here

    interpArg = cell(1,length(voltageVec));
    for ii = 1:length(voltageVec)
        interpArg{ii} = voltageVec(ii);
    end
    
    % 2DEG
    if nargin == 2
        interpArg = [interpArg, xx];
    % 2D
    elseif nargin == 3
        interpArg = [interpArg, zz, xx];
    end
end

