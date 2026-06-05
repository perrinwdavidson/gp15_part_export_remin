%%  calculate all stations
%   loop through all stations ::
for iStat = 1 : 1 : NUMSTAT

	%   get station number ::
	statNo = gp15_stations.stationNo(iStat);
	idxStat = gp15_flux.stationNo == statNo;

	%   get depth ::
	x = gp15_flux.depth(idxStat);

	%   get data ::
	yth = gp15_flux.th234(idxStat);
	yu = gp15_flux.u238(idxStat);
	for dataProduct = {'ecco', 'cglo', 'foam', 'glor', 'oras', 'grep', 'meanFull', 'meanNoFOAM'}
		dataProduct = dataProduct{1};
		% E4: eval removed; use dynamic struct field access
		% eval(['yw_' dataProduct ' = gp15_flux.w_' dataProduct '(idxStat);']);
		yw.(dataProduct) = gp15_flux.(['w_' dataProduct])(idxStat);
	end
	ygrad = gp15_flux.vertGrad(idxStat);

	%   get error ::
	nth = gp15_flux.uncertTh234(idxStat);  % assumed to be standard deviation
	nu = gp15_flux.uncertU238(idxStat);  % assumed to be standard deviation
	for dataProduct = {'ecco', 'cglo', 'foam', 'glor', 'oras', 'grep', 'meanFull', 'meanNoFOAM'}
		dataProduct = dataProduct{1};
		% E4: eval removed; use dynamic struct field access
		% eval(['nw_' dataProduct ' = gp15_flux.wErr_' dataProduct '(idxStat);']);  % this is the standard error
		nw.(dataProduct) = gp15_flux.(['wErr_' dataProduct])(idxStat);  % this is the standard error
	end
	ngrad = gp15_flux.vertGradError(idxStat);  % this is the noise term of the spline assumed to be the standard deviation

	%   calculate number of depths ::
	n = length(x);

	%   calculate dz ::
	% C2: inline dz computation extracted to calcLayerThickness(x)
	% dz1 = x;
	% dz0 = [0; dz1(1 : (end - 1))];
	% dz2 = [dz1(2 : end); dz1(end)];
	% bc(1, 1) = 2;
	% bc(2 : length(x), 1) = 1;
	% dz = ((dz2 - dz1) ./ 2) + (((dz1 - dz0) ./ 2) .* bc);
	dz = calcLayerThickness(x);

	%   calculate simple (1d) error ::
	%%% calculate ::
	error1d = cumsum((((LAMBDA * dz) .^ 2) .* ((nu .^ 2) + (nth .^ 2))), 'omitnan');  % this is the variance

	%%% calculate poc error ::
	% pocRatioUncert = gp15_flux.uncertPocTh234RatioRegress(idxStat) ./ (gp15_flux.pocTh234RatioRegress(idxStat) .^ 2);
	pocRatioUncert = (gp15_flux.uncertPocTh234RatioRegress(idxStat) .^ 2) ./ (gp15_flux.pocTh234RatioRegress(idxStat) .^ 2);
	pocRatioUncert(isnan(pocRatioUncert)) = 0;
	error1dPoc = (gp15_flux.pocFluxCumul1d(idxStat) .^ 2) .* ((error1d ./ (gp15_flux.th234FluxCumul1d(idxStat) .^ 2)) + pocRatioUncert);

	%%% calculate pon error ::
	% pnRatioUncert = gp15_flux.uncertPnTh234RatioRegress(idxStat) ./ (gp15_flux.pnTh234RatioRegress(idxStat) .^ 2);
	pnRatioUncert = (gp15_flux.uncertPnTh234RatioRegress(idxStat) .^ 2) ./ (gp15_flux.pnTh234RatioRegress(idxStat) .^ 2);
	pnRatioUncert(isnan(pnRatioUncert)) = 0;
	error1dPn = (gp15_flux.pnFluxCumul1d(idxStat) .^ 2) .* ((error1d ./ (gp15_flux.th234FluxCumul1d(idxStat) .^ 2)) + pnRatioUncert);

	%%% store ::
	gp15_flux.uncert1d(idxStat) = sqrt(error1d);  % this is the standard deviation now
	gp15_flux.uncert1dPoc(idxStat) = sqrt(error1dPoc);  % this is the standard deviation now
	gp15_flux.uncert1dPn(idxStat) = sqrt(error1dPn);  % this is the standard deviation now

	for dataProduct = {'ecco', 'cglo', 'foam', 'glor', 'oras', 'grep', 'meanFull', 'meanNoFOAM'}

		% get data product ::
		dataProduct = dataProduct{1};

		% calculate upwelling error ::
		% E4: eval removed; use dynamic struct field access
		% eval(['ywdp = yw_' dataProduct ';']);
		% eval(['nwdp = nw_' dataProduct ';']);
		ywdp = yw.(dataProduct);
		nwdp = nw.(dataProduct);
		% E1: relative-error form produces 0/0 = NaN when w or grad is zero; switch to absolute-error form
		% upwellError = ((ywdp .* ygrad) .^ 2) .* (((nwdp ./ ywdp) .^ 2) + ((ngrad ./ ygrad) .^ 2));
		upwellError = (ygrad .* dz) .^ 2 .* nwdp .^ 2 + (ywdp .* dz) .^ 2 .* ngrad .^ 2;
		% E2: NaN→0 before cumsum hides genuinely missing variance; propagate NaN instead
		% upwellError(isnan(upwellError)) = 0;
		% E2: 'omitnan' removed so NaN in upwellError propagates through cumulative sum
		% error2d = cumsum((((LAMBDA * dz) .^ 2) .* ((nu .^ 2) + (nth .^ 2)) + upwellError), 'omitnan');  % this is again variance
		error2d = cumsum((((LAMBDA * dz) .^ 2) .* ((nu .^ 2) + (nth .^ 2)) + upwellError));  % this is again variance

		% get th234 flux ::
		% E4: eval removed; use dynamic struct field access
		% eval(['flux2d = gp15_flux.th234FluxCumul2d_' dataProduct '(idxStat);']);
		flux2d = gp15_flux.(['th234FluxCumul2d_' dataProduct])(idxStat);
		% fluxError = error2d ./ flux2d;  % bug #5: divides by F_Th (not F_Th^2); gives Var/F instead of Var/F^2; inconsistent with 1D formulation
		fluxError = error2d ./ (flux2d .^ 2);
		fluxError(isnan(fluxError)) = 0;

		% get particulate flux ::
		% E4: eval removed; use dynamic struct field access
		% eval(['pocFlux = gp15_flux.pocFluxCumul2d_' dataProduct '(idxStat);']);
		% eval(['pnFlux = gp15_flux.pnFluxCumul2d_' dataProduct '(idxStat);']);
		pocFlux = gp15_flux.(['pocFluxCumul2d_' dataProduct])(idxStat);
		pnFlux = gp15_flux.(['pnFluxCumul2d_' dataProduct])(idxStat);

		% calculate particulate error ::
		% error2dPoc = (pocFlux .^ 2) .* ((fluxError .^ 2) + (pocRatioUncert .^ 2));
		% error2dPoc = (pocFlux .^ 2) .* ((fluxError .^ 2) + (pocRatioUncert));  % bug #5: fluxError is now Var/F^2; squaring gives Var^2/F^4 instead of Var/F^2
		error2dPoc = (pocFlux .^ 2) .* (fluxError + pocRatioUncert);
		% error2dPn = (pnFlux .^ 2) .* ((fluxError .^ 2) + (pnRatioUncert .^ 2));
		% error2dPn = (pnFlux .^ 2) .* ((fluxError .^ 2) + (pnRatioUncert));  % bug #5: same as above for PN
		error2dPn = (pnFlux .^ 2) .* (fluxError + pnRatioUncert);
		% E2: NaN→0 hides missing ratio/upwelling data; propagate NaN instead
		% error2dPoc(isnan(error2dPoc)) = 0;  % as all poc or pn = nan is ratio = nan
		% error2dPn(isnan(error2dPn)) = 0;

		% store ::
		% E4: eval removed; use dynamic struct field access
		% eval(['gp15_flux.uncertUpwell_' dataProduct '(idxStat) = sqrt(error2d);']);
		% eval(['gp15_flux.uncertUpwellCorrect_' dataProduct '(idxStat) = gp15_flux.uncertUpwell_' dataProduct '(idxStat) - gp15_flux.uncert1d(idxStat);']);
		% eval(['gp15_flux.uncertUpwellPoc_' dataProduct '(idxStat) = sqrt(error2dPoc);']);
		% eval(['gp15_flux.uncertUpwellPocCorrect_' dataProduct '(idxStat) = gp15_flux.uncertUpwellPoc_' dataProduct '(idxStat) - gp15_flux.uncert1dPoc(idxStat);']);
		% eval(['gp15_flux.uncertUpwellPn_' dataProduct '(idxStat) = sqrt(error2dPn);']);
		% eval(['gp15_flux.uncertUpwellPnCorrect_' dataProduct '(idxStat) = gp15_flux.uncertUpwellPn_' dataProduct '(idxStat) - gp15_flux.uncert1dPn(idxStat);']);
		gp15_flux.(['uncertUpwell_' dataProduct])(idxStat) = sqrt(error2d);
		% E5: sigma_2D - sigma_1D understates the correction uncertainty; correct derivation:
		% F_2D = F_1D + dF_upwell; errors independent -> Var[dF] = Var[F_2D] - Var[F_1D]
		% -> sigma_correction = sqrt(max(0, sigma_2D^2 - sigma_1D^2))
		% old: gp15_flux.(['uncertUpwellCorrect_' dataProduct])(idxStat) = gp15_flux.(['uncertUpwell_' dataProduct])(idxStat) - gp15_flux.uncert1d(idxStat);
		gp15_flux.(['uncertUpwellCorrect_' dataProduct])(idxStat) = sqrt(max(0, gp15_flux.(['uncertUpwell_' dataProduct])(idxStat) .^ 2 - gp15_flux.uncert1d(idxStat) .^ 2));
		gp15_flux.(['uncertUpwellPoc_' dataProduct])(idxStat) = sqrt(error2dPoc);
		% old: gp15_flux.(['uncertUpwellPocCorrect_' dataProduct])(idxStat) = gp15_flux.(['uncertUpwellPoc_' dataProduct])(idxStat) - gp15_flux.uncert1dPoc(idxStat);
		gp15_flux.(['uncertUpwellPocCorrect_' dataProduct])(idxStat) = sqrt(max(0, gp15_flux.(['uncertUpwellPoc_' dataProduct])(idxStat) .^ 2 - gp15_flux.uncert1dPoc(idxStat) .^ 2));
		gp15_flux.(['uncertUpwellPn_' dataProduct])(idxStat) = sqrt(error2dPn);
		% old: gp15_flux.(['uncertUpwellPnCorrect_' dataProduct])(idxStat) = gp15_flux.(['uncertUpwellPn_' dataProduct])(idxStat) - gp15_flux.uncert1dPn(idxStat);
		gp15_flux.(['uncertUpwellPnCorrect_' dataProduct])(idxStat) = sqrt(max(0, gp15_flux.(['uncertUpwellPn_' dataProduct])(idxStat) .^ 2 - gp15_flux.uncert1dPn(idxStat) .^ 2));

	end

	%   clear ::
	% C2: dz0/dz1/dz2/bc removed (now encapsulated in calcLayerThickness)
	% clear('dz', 'dz0', 'dz1', 'dz2', 'bc', 'fluxErrorStat');
	clear('dz', 'yw', 'nw', 'fluxErrorStat');

end

%%  calculate statistics of upwelling ensemble flux
%   calculate standard error of models ::
%%% th234 ::
% th234FluxData = [gp15_flux.th234FluxCumul2d_ecco, gp15_flux.th234FluxCumul2d_cglo, gp15_flux.th234FluxCumul2d_glor, gp15_flux.th234FluxCumul2d_oras];  % foam excluded (previously noted as anomalous near equator)
% th234FluxData = [gp15_flux.th234FluxCumul2d_ecco, gp15_flux.th234FluxCumul2d_cglo, gp15_flux.th234FluxCumul2d_foam, gp15_flux.th234FluxCumul2d_glor, gp15_flux.th234FluxCumul2d_oras, gp15_flux.th234FluxCumul2d_grep];  % GREP excluded: it is the arithmetic mean of CGLO/FOAM/GLOR/ORAS and is not an independent member
% modelErrorTh234Flux = std(th234FluxData, 1, 2) / sqrt(size(th234FluxData, 2));  % /sqrt(n) removed: models are not i.i.d. samples; inter-model std is the structural uncertainty
th234FluxData = [gp15_flux.th234FluxCumul2d_ecco, gp15_flux.th234FluxCumul2d_cglo, gp15_flux.th234FluxCumul2d_foam, gp15_flux.th234FluxCumul2d_glor, gp15_flux.th234FluxCumul2d_oras];
modelErrorTh234Flux = std(th234FluxData, 1, 2);

%%% poc ::
% pocFluxData = [gp15_flux.pocFluxCumul2d_ecco, gp15_flux.pocFluxCumul2d_cglo, gp15_flux.pocFluxCumul2d_glor, gp15_flux.pocFluxCumul2d_oras];  % foam excluded
% pocFluxData = [gp15_flux.pocFluxCumul2d_ecco, gp15_flux.pocFluxCumul2d_cglo, gp15_flux.pocFluxCumul2d_foam, gp15_flux.pocFluxCumul2d_glor, gp15_flux.pocFluxCumul2d_oras, gp15_flux.pocFluxCumul2d_grep];  % GREP excluded; /sqrt(n) removed
pocFluxData = [gp15_flux.pocFluxCumul2d_ecco, gp15_flux.pocFluxCumul2d_cglo, gp15_flux.pocFluxCumul2d_foam, gp15_flux.pocFluxCumul2d_glor, gp15_flux.pocFluxCumul2d_oras];
modelErrorPocFlux = std(pocFluxData, 1, 2);

%%% pn ::
% pnFluxData = [gp15_flux.pnFluxCumul2d_ecco, gp15_flux.pnFluxCumul2d_cglo, gp15_flux.pnFluxCumul2d_glor, gp15_flux.pnFluxCumul2d_oras];  % foam excluded
% pnFluxData = [gp15_flux.pnFluxCumul2d_ecco, gp15_flux.pnFluxCumul2d_cglo, gp15_flux.pnFluxCumul2d_foam, gp15_flux.pnFluxCumul2d_glor, gp15_flux.pnFluxCumul2d_oras, gp15_flux.pnFluxCumul2d_grep];  % GREP excluded; /sqrt(n) removed
pnFluxData = [gp15_flux.pnFluxCumul2d_ecco, gp15_flux.pnFluxCumul2d_cglo, gp15_flux.pnFluxCumul2d_foam, gp15_flux.pnFluxCumul2d_glor, gp15_flux.pnFluxCumul2d_oras];
modelErrorPnFlux = std(pnFluxData, 1, 2);

%   calculate total error, full quadrature combination of measurement and model-structural uncertainty ::
% gp15_flux.totalTh234FluxError = gp15_flux.uncert1d + sqrt((modelErrorTh234Flux .^ 2) + (gp15_flux.uncertUpwellCorrect_ecco .^ 2));  % original: linear addition assumes 100% correlation; correction term subtracts SDs not variances
gp15_flux.totalTh234FluxError = sqrt((gp15_flux.uncertUpwell_ecco .^ 2) + (modelErrorTh234Flux .^ 2));
% gp15_flux.totalPocFluxError = gp15_flux.uncert1dPoc + sqrt((modelErrorPocFlux .^ 2) + (gp15_flux.uncertUpwellPocCorrect_ecco .^ 2));  % original: same issue
gp15_flux.totalPocFluxError = sqrt((gp15_flux.uncertUpwellPoc_ecco .^ 2) + (modelErrorPocFlux .^ 2));
% gp15_flux.totalPnFluxError = gp15_flux.uncert1dPn + sqrt((modelErrorPnFlux .^ 2) + (gp15_flux.uncertUpwellPnCorrect_ecco .^ 2));  % original: same issue
gp15_flux.totalPnFluxError = sqrt((gp15_flux.uncertUpwellPn_ecco .^ 2) + (modelErrorPnFlux .^ 2));

%%  add named reporting-product total errors
%   ecco alias: σ_2d,ECCO + σ_model,full — backward-compat copy of the unsuffixed columns ::
gp15_flux.totalTh234FluxError_ecco = gp15_flux.totalTh234FluxError;
gp15_flux.totalPocFluxError_ecco   = gp15_flux.totalPocFluxError;
gp15_flux.totalPnFluxError_ecco    = gp15_flux.totalPnFluxError;

%   no-FOAM ensemble spread (n=4: ECCO, CGLO, GLOR, ORAS) ::
th234FluxData_noFOAM   = [gp15_flux.th234FluxCumul2d_ecco, gp15_flux.th234FluxCumul2d_cglo, gp15_flux.th234FluxCumul2d_glor, gp15_flux.th234FluxCumul2d_oras];
modelErrorTh234Flux_noFOAM = std(th234FluxData_noFOAM, 1, 2);
pocFluxData_noFOAM     = [gp15_flux.pocFluxCumul2d_ecco, gp15_flux.pocFluxCumul2d_cglo, gp15_flux.pocFluxCumul2d_glor, gp15_flux.pocFluxCumul2d_oras];
modelErrorPocFlux_noFOAM   = std(pocFluxData_noFOAM, 1, 2);
pnFluxData_noFOAM      = [gp15_flux.pnFluxCumul2d_ecco, gp15_flux.pnFluxCumul2d_cglo, gp15_flux.pnFluxCumul2d_glor, gp15_flux.pnFluxCumul2d_oras];
modelErrorPnFlux_noFOAM    = std(pnFluxData_noFOAM, 1, 2);

%   meanFull: σ_2d,meanFull + σ_model,full (n=5 ensemble) ::
gp15_flux.totalTh234FluxError_meanFull = sqrt((gp15_flux.uncertUpwell_meanFull .^ 2)    + (modelErrorTh234Flux .^ 2));
gp15_flux.totalPocFluxError_meanFull   = sqrt((gp15_flux.uncertUpwellPoc_meanFull .^ 2) + (modelErrorPocFlux .^ 2));
gp15_flux.totalPnFluxError_meanFull    = sqrt((gp15_flux.uncertUpwellPn_meanFull .^ 2)  + (modelErrorPnFlux .^ 2));

%   meanNoFOAM: σ_2d,meanNoFOAM + σ_model,noFOAM (n=4 ensemble) ::
gp15_flux.totalTh234FluxError_meanNoFOAM = sqrt((gp15_flux.uncertUpwell_meanNoFOAM .^ 2)    + (modelErrorTh234Flux_noFOAM .^ 2));
gp15_flux.totalPocFluxError_meanNoFOAM   = sqrt((gp15_flux.uncertUpwellPoc_meanNoFOAM .^ 2) + (modelErrorPocFlux_noFOAM .^ 2));
gp15_flux.totalPnFluxError_meanNoFOAM    = sqrt((gp15_flux.uncertUpwellPn_meanNoFOAM .^ 2)  + (modelErrorPnFlux_noFOAM .^ 2));

%%  end subroutine
