function dz = calcLayerThickness(z)
%%  half-open midpoint-rule cell widths
%   surface cell (i=1) extends to z=0 (bc=2); all others use centered half-spacings.
dz1 = z;
dz0 = [0; dz1(1 : end - 1)];
dz2 = [dz1(2 : end); dz1(end)];
bc = ones(size(z));
bc(1) = 2;
dz = ((dz2 - dz1) ./ 2) + (((dz1 - dz0) ./ 2) .* bc);
end
