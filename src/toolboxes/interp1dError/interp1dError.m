function [vq, vqStd] = interp1dError(x, v, vStd, xq0)

	% check data input ::
	blank = NaN(size(xq0));
	if (sum(isnan(x)) > 0) || (sum(isnan(v)) > 0)
		vq = blank; 
		vqStd = blank; 
		return
	end

	% bin data ::
	[~, k0] = histc(xq0, x); 
	n = length(x);
	k0(k0 == n) = n - 1;

	% handle domain ::
	xq = xq0(k0 ~= 0);
	k = k0(k0 ~= 0); 
	
	% make interpolation polynomials ::
	p1 = (xq - x(k+1)) ./ (x(k) - x(k+1));
	p2 = (xq - x(k)) ./ (x(k+1) - x(k));

	% interpolate ::
	vq0 = (v(k) .* p1) + (v(k+1) .* p2);
	vqStd0 = sqrt(((vStd(k) .* p1) .^ 2) + ((vStd(k+1) .* p2) .^ 2));

	% pad where no extrapolation ::
	vq = blank;
	vqStd = blank;
	vq(k0 ~= 0) = vq0;
	vqStd(k0 ~= 0) = vqStd0;

end
