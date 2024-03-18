%%  calculate poc and pn flux
%   1d ::
gp15_flux.pocFluxCumul1d = gp15_flux.th234FluxCumul1d .* gp15_flux.pocTh234RatioRegress;
gp15_flux.pnFluxCumul1d = gp15_flux.th234FluxCumul1d .* gp15_flux.pnTh234RatioRegress;

%  2d ::
for dataProduct = {'ecco', 'cglo', 'foam', 'glor', 'oras'}

	% get data product ::
	dataProduct = dataProduct{1}; 

	% get fluxes ::
	eval(['gp15_flux.pocFluxCumul2d_' dataProduct ' = gp15_flux.th234FluxCumul2d_' dataProduct ' .* gp15_flux.pocTh234RatioRegress;']);
	eval(['gp15_flux.pnFluxCumul2d_' dataProduct ' = gp15_flux.th234FluxCumul2d_' dataProduct ' .* gp15_flux.pnTh234RatioRegress;']);

end

%%  end subroutine
