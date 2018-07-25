function interpArg = getInterpolantArgument( voltageVec, xx, zz )
%GETINTERPOLANTARGUMENT Summary of this function goes here
%   Detailed explanation goes here

    interpArg = cell(1,length(voltageVec)+2);
    for ii = 1:length(voltageVec)
        interpArg{ii} = voltageVec(ii);
    end
    
    % 2DEG
    if nargin == 2
        interpArg{length(voltageVec)+1} = xx;
        interpArg(length(voltageVec)+2) = [];
    % 2D
    elseif nargin == 3
        interpArg{length(voltageVec)+1} = zz;
        interpArg{length(voltageVec)+2} = xx;
    end
end

