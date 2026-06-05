%%  make depth arrays
%   set tolerance ::
tol = 1e-6;

%   100m ::
%   bug: old findDepthIndex used min() which returned a single scalar index;
%        gp15_flux has one row per (station, depth) pair so depth 100 appears
%        NUMSTAT times — the single-index approach silently returned only one row.
% gp15_flux100 = gp15_flux(gp15_flux.depth == 100, :);
% [~, idx100, ~] = unique(gp15_flux100.latitude);
% gp15_flux100 = gp15_flux100(idx100, :);
% idx100 = findDepthIndex(gp15_flux.depth, 100, tol, '100');  % bug: returns 1 index, not NUMSTAT
idx100 = find(abs(gp15_flux.depth - 100) <= tol);
gp15_flux100 = gp15_flux(idx100, :);
assert(height(gp15_flux100) == NUMSTAT, ...
       'Expected %d rows at 100 m depth; found %d. Check that 100 m was inserted at interpolation.', ...
       NUMSTAT, height(gp15_flux100));

%   ppz ::
%   bug: findDepthIndex called with a NUMSTAT×1 vector for target_depth;
%        depth_array - target_depth is a dimension mismatch (N×1 vs NUMSTAT×1)
%        and min() would return a single index even if it did not error.
% [~, idxPpz0] = intersect(gp15_flux.depth, gp15_stations.depthPpz);
% gp15_fluxPpz = gp15_flux(idxPpz0, :);
% [~, idxPpz, ~] = unique(gp15_fluxPpz.latitude);
% gp15_fluxPpz = gp15_fluxPpz(idxPpz, :);
% idxPpz = findDepthIndex(gp15_flux.depth, gp15_stations.depthPpz, tol, 'PPZ');  % bug: vector target fails
idxPpz = false(height(gp15_flux), 1);
for iPpz = 1 : 1 : NUMSTAT
    snPpz  = gp15_stations.stationNo(iPpz);
    ppzDep = gp15_stations.depthPpz(iPpz);
    idxPpz = idxPpz | ((gp15_flux.stationNo == snPpz) & (abs(gp15_flux.depth - ppzDep) <= tol));
end
gp15_fluxPpz = gp15_flux(idxPpz, :);
assert(height(gp15_fluxPpz) == NUMSTAT, ...
       'Expected %d rows at PPZ depth; found %d. Check PPZ depth insertion in stationInterp.', ...
       NUMSTAT, height(gp15_fluxPpz));

%   ppz + 100 ::
% [~, idx100Ppz0] = intersect(gp15_flux.depth, gp15_stations.depthPpz + 100);
% gp15_flux100Ppz = gp15_flux(idx100Ppz0, :);
% [~, idx100Ppz, ~] = unique(gp15_flux100Ppz.latitude);
% gp15_flux100Ppz = gp15_flux100Ppz(idx100Ppz, :);
% idx100Ppz = findDepthIndex(gp15_flux.depth, gp15_stations.depthPpz + 100, tol, 'PPZ+100');  % bug: vector target fails
idx100Ppz = false(height(gp15_flux), 1);
for iPpz = 1 : 1 : NUMSTAT
    snPpz  = gp15_stations.stationNo(iPpz);
    ppzDep = gp15_stations.depthPpz(iPpz) + 100;
    idx100Ppz = idx100Ppz | ((gp15_flux.stationNo == snPpz) & (abs(gp15_flux.depth - ppzDep) <= tol));
end
gp15_flux100Ppz = gp15_flux(idx100Ppz, :);
assert(height(gp15_flux100Ppz) == NUMSTAT, ...
       'Expected %d rows at PPZ+100 m depth; found %d. Check PPZ+100 depth insertion in stationInterp.', ...
       NUMSTAT, height(gp15_flux100Ppz));

%%  end subroutine
