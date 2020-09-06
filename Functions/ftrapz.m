function z = ftrapz(y, dim)

    nshifts = 0;

      dim = min(ndims(y)+1, dim);
      perm = [dim:max(ndims(y),dim) 1:dim-1];

    % Trapezoid sum computed with vector-matrix multiply.
    z = sum((y(1:end-1,:) + y(2:end,:)), 1)/2;

    siz = size(y);
    siz(1) = 1;
    z = reshape(z,[ones(1,nshifts),siz]);
    if ~isempty(perm) && ~isscalar(z)
        z = ipermute(z,perm);
    end
end

