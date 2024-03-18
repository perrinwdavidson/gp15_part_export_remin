%%  do model calculations
%   get station number ::
statNo = gp15_stations.stationNo(iStat);

%   specify station data array ::
statProfile = gp15_inputs(gp15_inputs.stationNo == statNo, :);

%   calculate dz ::
dz1 = statProfile.depth;
dz0 = [0; dz1(1 : (end - 1))];
dz2 = [dz1(2 : end); dz1(end)];
bc(1, 1) = 2;
bc(2 : length(statProfile.depth), 1) = 1;
statProfile.dz = ((dz2 - dz1) ./ 2) + (((dz1 - dz0) ./ 2) .* bc);

%   calculate radionuclide inventories ::
statProfile.inventoryU238 = statProfile.u238 .* statProfile.dz;
statProfile.inventoryTh234 = statProfile.th234 .* statProfile.dz;

%   calculate th-234 deficit ::
statProfile.th234Deficit = statProfile.inventoryU238 - statProfile.inventoryTh234;

%   calculate upwelling correction ::
%%% loop through all model products ::
for dataProduct = {'ecco', 'cglo', 'foam', 'glor', 'oras'}

	%   get data product name ::
	dataProduct = dataProduct{1}; 

	%   get upwelling data ::
	eval(['statProfile.upwell_' dataProduct ' = statProfile.w_' dataProduct ' .* statProfile.vertGrad .* statProfile.dz;']);

	%   calculate th-234 flux per layer ::
	if strcmp(dataProduct, 'ecco')
		statProfile.th234FluxPerLayer1d = statProfile.th234Deficit .* LAMBDA;
	end
	eval(['statProfile.th234FluxPerLayer2d_' dataProduct ' = statProfile.th234Deficit .* LAMBDA + statProfile.upwell_' dataProduct ';']);

	%   calculate cumulative fluxes ::
	for iDepth = 1 : 1 : length(statProfile.depth)
		if iDepth == 1

			%   th234 flux ::
			if strcmp(dataProduct, 'ecco')
				statProfile.th234FluxCumul1d(iDepth) = statProfile.th234FluxPerLayer1d(iDepth);
			end
			eval(['statProfile.th234FluxCumul2d_' dataProduct '(iDepth) = statProfile.th234FluxPerLayer2d_' dataProduct '(iDepth);']);

		else

			%   th234 flux ::
			if strcmp(dataProduct, 'ecco')
				statProfile.th234FluxCumul1d(iDepth) = statProfile.th234FluxPerLayer1d(iDepth) + statProfile.th234FluxCumul1d(iDepth-1);
			end
			eval(['statProfile.th234FluxCumul2d_' dataProduct '(iDepth) = statProfile.th234FluxPerLayer2d_' dataProduct '(iDepth) + statProfile.th234FluxCumul2d_' dataProduct '(iDepth-1);']);

		end
	end

end

%   store data ::
if iStat == 1

	gp15_flux = statProfile;

else

	gp15_flux = vertcat(gp15_flux, statProfile);

end

%   clear data ::
clear('dz0', 'dz1', 'dz2', 'bc'); 

%%  end subroutine
