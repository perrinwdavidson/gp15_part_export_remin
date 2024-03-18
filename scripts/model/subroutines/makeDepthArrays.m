%%  make depth arrays
%   100m ::
gp15_flux100 = gp15_flux(gp15_flux.depth == 100, :); 
[~, idx100, ~] = unique(gp15_flux100.latitude);
gp15_flux100 = gp15_flux100(idx100, :); 

%   ppz ::
[~, idxPpz0] = intersect(gp15_flux.depth, gp15_stations.depthPpz); 
gp15_fluxPpz = gp15_flux(idxPpz0, :); 
[~, idxPpz, ~] = unique(gp15_fluxPpz.latitude);
gp15_fluxPpz = gp15_fluxPpz(idxPpz, :); 

%   ppz + 100 ::
[~, idx100Ppz0] = intersect(gp15_flux.depth, gp15_stations.depthPpz + 100); 
gp15_flux100Ppz = gp15_flux(idx100Ppz0, :); 
[~, idx100Ppz, ~] = unique(gp15_flux100Ppz.latitude);
gp15_flux100Ppz = gp15_flux100Ppz(idx100Ppz, :); 

%%  end subroutine
