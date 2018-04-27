function [ gMG ] = convertNOtoMG( vNO, nx, ny )
%CONVERTNOTOMG converts a vector in natural ordering (NO) to a meshgrid
%(MG) format.
%   nx and ny give the number of columns and rows respectively of the
%   desired meshgrid format.  While we don't explicitly check for it.. It
%   is implied that nx*ny = length(wfNO)

    gMG = reshape(vNO,[nx,ny]).';
end

