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
	for dataProduct = {'ecco', 'cglo', 'foam', 'glor', 'oras'}
		dataProduct = dataProduct{1}; 
		eval(['yw_' dataProduct ' = gp15_flux.w_' dataProduct '(idxStat);']); 
	end
	ygrad = gp15_flux.vertGrad(idxStat); 

	%   get error ::
	nth = gp15_flux.uncertTh234(idxStat);  % assumed to be standard deviation 
	nu = gp15_flux.uncertU238(idxStat);  % assumed to be standard deviation 
	for dataProduct = {'ecco', 'cglo', 'foam', 'glor', 'oras'}
		dataProduct = dataProduct{1}; 
		eval(['nw_' dataProduct ' = gp15_flux.wErr_' dataProduct '(idxStat);']);  % this is the standard error
	end
	ngrad = gp15_flux.vertGradError(idxStat);  % this is the noise term of the spline assumed to be the standard deviation 

	%   calculate number of depths ::
	n = length(x); 

	%   calculate dz ::
	dz1 = x;
	dz0 = [0; dz1(1 : (end - 1))];
	dz2 = [dz1(2 : end); dz1(end)];
	bc(1, 1) = 2;
	bc(2 : length(x), 1) = 1;
	dz = ((dz2 - dz1) ./ 2) + (((dz1 - dz0) ./ 2) .* bc);

	%   calculate simple (1d) error ::
	%%% calculate ::
	error1d = cumsum((((LAMBDA * dz) .^ 2) .* ((nu .^ 2) + (nth .^ 2))), 'omitnan');  % this is the variance

	%%% calculate poc error ::
	pocRatioUncert = gp15_flux.uncertPocTh234RatioRegress(idxStat) ./ (gp15_flux.pocTh234RatioRegress(idxStat) .^ 2);
	pocRatioUncert(isnan(pocRatioUncert)) = 0; 
	error1dPoc = (gp15_flux.pocFluxCumul1d(idxStat) .^ 2) .* ((error1d ./ (gp15_flux.th234FluxCumul1d(idxStat) .^ 2)) + pocRatioUncert); 
	
	%%% calculate poc error ::
	pnRatioUncert = gp15_flux.uncertPnTh234RatioRegress(idxStat) ./ (gp15_flux.pnTh234RatioRegress(idxStat) .^ 2);
	pnRatioUncert(isnan(pnRatioUncert)) = 0; 
	error1dPn = (gp15_flux.pnFluxCumul1d(idxStat) .^ 2) .* ((error1d ./ (gp15_flux.th234FluxCumul1d(idxStat) .^ 2)) + pnRatioUncert); 

	%%% store ::
	gp15_flux.uncert1d(idxStat) = sqrt(error1d);  % this is the standard deviation now
	gp15_flux.uncert1dPoc(idxStat) = sqrt(error1dPoc);  % this is the standard deviation now
	gp15_flux.uncert1dPn(idxStat) = sqrt(error1dPn);  % this is the standard deviation now

	for dataProduct = {'ecco', 'cglo', 'foam', 'glor', 'oras'}

		% get data product ::
		dataProduct = dataProduct{1}; 

		% calculate upwelling error ::
		eval(['ywdp = yw_' dataProduct ';']); 
		upwellError = ((ywdp .* ygrad) .* (((nw_ecco ./ ywdp) .^ 2) + ((ngrad ./ ygrad) .^ 2)));
		upwellError(isnan(upwellError)) = 0; 
		error2d = cumsum((((LAMBDA * dz) .^ 2) .* ((nu .^ 2) + (nth .^ 2)) + upwellError), 'omitnan');  % this is again variance

		% get th234 flux ::
		eval(['flux2d = gp15_flux.th234FluxCumul2d_' dataProduct '(idxStat);']); 
		fluxError = error2d ./ flux2d; 
		fluxError(isnan(fluxError)) = 0; 

		% get particulate flux ::
		eval(['pocFlux = gp15_flux.pocFluxCumul2d_' dataProduct '(idxStat);']);
		eval(['pnFlux = gp15_flux.pnFluxCumul2d_' dataProduct '(idxStat);']);

		% calculate particulate error ::
		error2dPoc = (pocFlux .^ 2) .* ((fluxError .^ 2) + (pocRatioUncert .^ 2)); 
		error2dPn = (pnFlux .^ 2) .* ((fluxError .^ 2) + (pnRatioUncert .^ 2)); 
		error2dPoc(isnan(error2dPoc)) = 0;  % as all poc or pn = nan is ratio = nan
		error2dPn(isnan(error2dPn)) = 0; 

		% store ::
		eval(['gp15_flux.uncertUpwell_' dataProduct '(idxStat) = sqrt(error2d);']);  
		eval(['gp15_flux.uncertUpwellCorrect_' dataProduct '(idxStat) = gp15_flux.uncertUpwell_' dataProduct '(idxStat) - gp15_flux.uncert1d(idxStat);']);  
		eval(['gp15_flux.uncertUpwellPoc_' dataProduct '(idxStat) = sqrt(error2dPoc);']);  
		eval(['gp15_flux.uncertUpwellPocCorrect_' dataProduct '(idxStat) = gp15_flux.uncertUpwellPoc_' dataProduct '(idxStat) - gp15_flux.uncert1dPoc(idxStat);']);  
		eval(['gp15_flux.uncertUpwellPn_' dataProduct '(idxStat) = sqrt(error2dPn);']);  
		eval(['gp15_flux.uncertUpwellPnCorrect_' dataProduct '(idxStat) = gp15_flux.uncertUpwellPn_' dataProduct '(idxStat) - gp15_flux.uncert1dPn(idxStat);']);  

	end

	%   clear ::
	clear('dz', 'dz0', 'dz1', 'dz2', 'bc', 'fluxErrorStat'); 
    
end

%%  calculate statistics of upwelling ensemble flux
%   calculate standard error of models ::
%   n.b.: foam is anamalous near equator so we neglect. 
%%% th234 ::
th234FluxData = [gp15_flux.th234FluxCumul2d_ecco, gp15_flux.th234FluxCumul2d_cglo, gp15_flux.th234FluxCumul2d_glor, gp15_flux.th234FluxCumul2d_oras];  % gp15_flux.th234FluxCumul2d_foam
modelErrorTh234Flux = std(th234FluxData, 1, 2) / sqrt(size(th234FluxData, 2));  % normalize by n (assume limit of N large)

%%% poc ::
pocFluxData = [gp15_flux.pocFluxCumul2d_ecco, gp15_flux.pocFluxCumul2d_cglo, gp15_flux.pocFluxCumul2d_glor, gp15_flux.pocFluxCumul2d_oras];  % gp15_flux.pocFluxCumul2d_foam
modelErrorPocFlux = std(pocFluxData, 1, 2) / sqrt(size(pocFluxData, 2)); 

%%% poc ::
pnFluxData = [gp15_flux.pnFluxCumul2d_ecco, gp15_flux.pnFluxCumul2d_cglo, gp15_flux.pnFluxCumul2d_glor, gp15_flux.pnFluxCumul2d_oras];  % gp15_flux.pnFluxCumul2d_foam
modelErrorPnFlux = std(pnFluxData, 1, 2) / sqrt(size(pnFluxData, 2)); 

%   calculate total error ::
gp15_flux.totalTh234FluxError = sqrt((modelErrorTh234Flux .^ 2) + (gp15_flux.uncertUpwellCorrect_ecco .^ 2)); 
gp15_flux.totalPocFluxError = sqrt((modelErrorPocFlux .^ 2) + (gp15_flux.uncertUpwellPocCorrect_ecco .^ 2)); 
gp15_flux.totalPnFluxError = sqrt((modelErrorPnFlux .^ 2) + (gp15_flux.uncertUpwellPnCorrect_ecco .^ 2)); 

%%  end subroutine
