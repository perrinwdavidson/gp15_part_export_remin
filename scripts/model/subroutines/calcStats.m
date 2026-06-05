%%  calculate stats
%   configurable reporting products ::
%   fluxProduct:   model product used for the reported 2D flux magnitude
%   uncertProduct: total-error column used for the reported uncertainty
%   valid values: 'ecco', 'cglo', 'foam', 'glor', 'oras', 'grep', 'meanFull', 'meanNoFOAM'
%   example: fluxProduct='ecco', uncertProduct='meanFull' reports ECCO magnitude with full-ensemble uncertainty
fluxProduct   = 'meanFull';
uncertProduct = 'meanFull';

%   all depth data of interest ::
%%% preallocate ::
depthData = []; 

%%% loop through stations ::
for iStat = 1 : 1 : NUMSTAT

	%   station number ::
	statNo = gp15_stations.stationNo(iStat); 
	
	%   zone ::
	statRegion = gp15_flux100.region(gp15_flux100.stationNo == statNo); 

	%   latitude ::
	statLat = round(gp15_stations.latitude(iStat), 1); 

	%   ppz depth ::
	statPpz = round(gp15_stations.depthPpz(iStat), 0); 

	%   100 m thorium flux ::
	%%% 1d ::
	stat234Th100mFlux = round(gp15_flux100.th234FluxCumul1d(gp15_flux100.stationNo == statNo), 2, 'significant'); 
	uncertStat234Th100mFlux = round(gp15_flux100.uncert1d(gp15_flux100.stationNo == statNo), 1, 'significant');

	%%% 2d ::
	% m4: old: hardcoded equatorial station list — silent failure if station numbers change
	% if isempty(intersect(statNo, [27, 29, 31]))  % <- THIS SHOULD BE MORE DEPENDENT ON STATISTICS
	if abs(statLat) > 5   % outside equatorial band (|lat| <= 5 deg)
		stat234Th100mFluxUpwell = NaN;
		uncertStat234Th100mFluxUpwell = NaN;
	else
		stat234Th100mFluxUpwell = round(gp15_flux100.(['th234FluxCumul2d_' fluxProduct])(gp15_flux100.stationNo == statNo), 2, 'significant');
		uncertStat234Th100mFluxUpwell = round(gp15_flux100.(['totalTh234FluxError_' uncertProduct])(gp15_flux100.stationNo == statNo), 1, 'significant');
	end

	%   ppz thorium flux ::
	%%% 1d ::
	stat234ThPpzFlux = round(gp15_fluxPpz.th234FluxCumul1d(gp15_fluxPpz.stationNo == statNo), 2, 'significant'); 
	uncertStat234ThPpzFlux = round(gp15_fluxPpz.uncert1d(gp15_fluxPpz.stationNo == statNo), 1, 'significant');

	%%% 2d ::
	% m4: old: hardcoded equatorial station list — silent failure if station numbers change
	% if isempty(intersect(statNo, [27, 29, 31]))
	if abs(statLat) > 5   % outside equatorial band (|lat| <= 5 deg)
		stat234ThPpzFluxUpwell = NaN;
		uncertStat234ThPpzFluxUpwell = NaN;
	else
		stat234ThPpzFluxUpwell = round(gp15_fluxPpz.(['th234FluxCumul2d_' fluxProduct])(gp15_fluxPpz.stationNo == statNo), 2, 'significant');
		uncertStat234ThPpzFluxUpwell = round(gp15_fluxPpz.(['totalTh234FluxError_' uncertProduct])(gp15_fluxPpz.stationNo == statNo), 1, 'significant');
	end
	
	%   100 m below ppz thorium flux ::
	%%% 1d ::
	stat234Th100PpzFlux = round(gp15_flux100Ppz.th234FluxCumul1d(gp15_flux100Ppz.stationNo == statNo), 2, 'significant'); 
	uncertStat234Th100PpzFlux = round(gp15_flux100Ppz.uncert1d(gp15_flux100Ppz.stationNo == statNo), 1, 'significant');

	%%% 2d ::
	% m4: old: hardcoded equatorial station list — silent failure if station numbers change
	% if isempty(intersect(statNo, [27, 29, 31]))
	if abs(statLat) > 5   % outside equatorial band (|lat| <= 5 deg)
		stat234Th100PpzFluxUpwell = NaN;
		uncertStat234Th100PpzFluxUpwell = NaN;
	else
		stat234Th100PpzFluxUpwell = round(gp15_flux100Ppz.(['th234FluxCumul2d_' fluxProduct])(gp15_flux100Ppz.stationNo == statNo), 2, 'significant');
		uncertStat234Th100PpzFluxUpwell = round(gp15_flux100Ppz.(['totalTh234FluxError_' uncertProduct])(gp15_flux100Ppz.stationNo == statNo), 1, 'significant');
	end

	%   100 m poc flux ::
	%%% 1d ::
	statPoc100mFlux = round(gp15_flux100.pocFluxCumul1d(gp15_flux100.stationNo == statNo), 2, 'significant'); 
	uncertStatPoc100mFlux = round(gp15_flux100.uncert1dPoc(gp15_flux100.stationNo == statNo), 1, 'significant');

	%%% 2d ::
	% m4: old: hardcoded equatorial station list — silent failure if station numbers change
	% if isempty(intersect(statNo, [27, 29, 31]))
	if abs(statLat) > 5   % outside equatorial band (|lat| <= 5 deg)
		statPoc100mFluxUpwell = NaN;
		uncertStatPoc100mFluxUpwell = NaN;
	else
		statPoc100mFluxUpwell = round(gp15_flux100.(['pocFluxCumul2d_' fluxProduct])(gp15_flux100.stationNo == statNo), 2, 'significant');
		uncertStatPoc100mFluxUpwell = round(gp15_flux100.(['totalPocFluxError_' uncertProduct])(gp15_flux100.stationNo == statNo), 1, 'significant');
	end

	%   ppz poc flux ::
	%%% 1d ::
	statPocPpzFlux = round(gp15_fluxPpz.pocFluxCumul1d(gp15_fluxPpz.stationNo == statNo), 2, 'significant'); 
	uncertStatPocPpzFlux = round(gp15_fluxPpz.uncert1dPoc(gp15_fluxPpz.stationNo == statNo), 1, 'significant');

	%%% 2d ::
	% m4: old: hardcoded equatorial station list — silent failure if station numbers change
	% if isempty(intersect(statNo, [27, 29, 31]))
	if abs(statLat) > 5   % outside equatorial band (|lat| <= 5 deg)
		statPocPpzFluxUpwell = NaN;
		uncertStatPocPpzFluxUpwell = NaN;
	else
		statPocPpzFluxUpwell = round(gp15_fluxPpz.(['pocFluxCumul2d_' fluxProduct])(gp15_fluxPpz.stationNo == statNo), 2, 'significant');
		uncertStatPocPpzFluxUpwell = round(gp15_fluxPpz.(['totalPocFluxError_' uncertProduct])(gp15_fluxPpz.stationNo == statNo), 1, 'significant');
	end

	%   100 m below ppz poc flux ::
	%%% 1d ::
	statPoc100PpzFlux = round(gp15_flux100Ppz.pocFluxCumul1d(gp15_flux100Ppz.stationNo == statNo), 2, 'significant'); 
	uncertStatPoc100PpzFlux = round(gp15_flux100Ppz.uncert1dPoc(gp15_flux100Ppz.stationNo == statNo), 1, 'significant');

	%%% 2d ::
	% m4: old: hardcoded equatorial station list — silent failure if station numbers change
	% if isempty(intersect(statNo, [27, 29, 31]))
	if abs(statLat) > 5   % outside equatorial band (|lat| <= 5 deg)
		statPoc100PpzFluxUpwell = NaN;
		uncertStatPoc100PpzFluxUpwell = NaN;
	else
		statPoc100PpzFluxUpwell = round(gp15_flux100Ppz.(['pocFluxCumul2d_' fluxProduct])(gp15_flux100Ppz.stationNo == statNo), 2, 'significant');
		uncertStatPoc100PpzFluxUpwell = round(gp15_flux100Ppz.(['totalPocFluxError_' uncertProduct])(gp15_flux100Ppz.stationNo == statNo), 1, 'significant');
	end

	%   100 m pn flux ::
	%%% 1d ::
	statPn100mFlux = round(gp15_flux100.pnFluxCumul1d(gp15_flux100.stationNo == statNo), 2, 'significant'); 
	uncertStatPn100mFlux = round(gp15_flux100.uncert1dPn(gp15_flux100.stationNo == statNo), 1, 'significant');

	%%% 2d ::
	% m4: old: hardcoded equatorial station list — silent failure if station numbers change
	% if isempty(intersect(statNo, [27, 29, 31]))
	if abs(statLat) > 5   % outside equatorial band (|lat| <= 5 deg)
		statPn100mFluxUpwell = NaN;
		uncertStatPn100mFluxUpwell = NaN;
	else
		statPn100mFluxUpwell = round(gp15_flux100.(['pnFluxCumul2d_' fluxProduct])(gp15_flux100.stationNo == statNo), 2, 'significant');
		uncertStatPn100mFluxUpwell = round(gp15_flux100.(['totalPnFluxError_' uncertProduct])(gp15_flux100.stationNo == statNo), 1, 'significant');
	end

	%   ppz pn flux ::
	%%% 1d ::
	statPnPpzFlux = round(gp15_fluxPpz.pnFluxCumul1d(gp15_fluxPpz.stationNo == statNo), 2, 'significant'); 
	uncertStatPnPpzFlux = round(gp15_fluxPpz.uncert1dPn(gp15_fluxPpz.stationNo == statNo), 1, 'significant');

	%%% 2d ::
	% m4: old: hardcoded equatorial station list — silent failure if station numbers change
	% if isempty(intersect(statNo, [27, 29, 31]))
	if abs(statLat) > 5   % outside equatorial band (|lat| <= 5 deg)
		statPnPpzFluxUpwell = NaN;
		uncertStatPnPpzFluxUpwell = NaN;
	else
		statPnPpzFluxUpwell = round(gp15_fluxPpz.(['pnFluxCumul2d_' fluxProduct])(gp15_fluxPpz.stationNo == statNo), 2, 'significant');
		uncertStatPnPpzFluxUpwell = round(gp15_fluxPpz.(['totalPnFluxError_' uncertProduct])(gp15_fluxPpz.stationNo == statNo), 1, 'significant');
	end

	%   make array ::
	depthData = [depthData; [statRegion, statNo, statLat, statPpz, ...
		   	         stat234Th100mFlux, uncertStat234Th100mFlux, stat234Th100mFluxUpwell, uncertStat234Th100mFluxUpwell, ...
			         stat234ThPpzFlux, uncertStat234ThPpzFlux, stat234ThPpzFluxUpwell, uncertStat234ThPpzFluxUpwell, ...
			         stat234Th100PpzFlux, uncertStat234Th100PpzFlux, stat234Th100PpzFluxUpwell, uncertStat234Th100PpzFluxUpwell, ...
			         statPoc100mFlux, uncertStatPoc100mFlux, statPoc100mFluxUpwell, uncertStatPoc100mFluxUpwell, ...
			         statPocPpzFlux, uncertStatPocPpzFlux, statPocPpzFluxUpwell, uncertStatPocPpzFluxUpwell, ...
			         statPoc100PpzFlux, uncertStatPoc100PpzFlux, statPoc100PpzFluxUpwell, uncertStatPoc100PpzFluxUpwell, ...
			         statPn100mFlux, uncertStatPn100mFlux, statPn100mFluxUpwell, uncertStatPn100mFluxUpwell, ...
			         statPnPpzFlux, uncertStatPnPpzFlux, statPnPpzFluxUpwell, uncertStatPnPpzFluxUpwell]]; 

end

%%% finalize table ::
depthData = array2table(depthData); 
depthData.Properties.VariableNames = {'Region', 'Station', 'Latitude', 'PPZ Depth (m)', ...
			              'P_Th(100) (dpm m-2 d-1)', 'P_Th_Std(100) (dpm m-2 d-1)', ...
				      'Upwelling-adjusted P_Th(100) (dpm m-2 d-1)', 'Upwelling-adjusted P_Th_Std(100) (dpm m-2 d-1)', ...
				      'P_Th(E_z) (dpm m-2 d-1)', 'P_Th_Std(E_z) (dpm m-2 d-1)', ...
			   	      'Upwelling-adjusted P_Th(E_z) (dpm m-2 d-1)', 'Upwelling-adjusted P_Th_Std(E_z) (dpm m-2 d-1)', ...
				      'P_Th(E_z+100) (dpm m-2 d-1)', 'P_Th_Std(E_z+100) (dpm m-2 d-1)', ...
			   	      'Upwelling-adjusted P_Th(E_z+100) (dpm m-2 d-1)', 'Upwelling-adjusted P_Th_Std(E_z+100) (dpm m-2 d-1)', ...
				      'P_POC(100) (mmol C m-2 d-1)', 'P_POC_Std(100) (mmol C m-2 d-1)', ...
				      'Upwelling-adjusted P_POC(100) (mmol C m-2 d-1)', 'Upwelling-adjusted P_POC_Std(100) (mmol C m-2 d-1)', ...
				      'P_POC(E_z) (mmol C m-2 d-1)', 'P_POC_Std(E_z) (mmol C m-2 d-1)', ...
			   	      'Upwelling-adjusted P_POC(E_z) (mmol C m-2 d-1)', 'Upwelling-adjusted P_POC_Std(E_z) (mmol C m-2 d-1)', ...
				      'P_POC(E_z+100) (mmol C m-2 d-1)', 'P_POC_Std(E_z+100) (mmol C m-2 d-1)', ...
			   	      'Upwelling-adjusted P_POC(E_z+100) (mmol C m-2 d-1)', 'Upwelling-adjusted P_POC_Std(E_z+100) (mmol C m-2 d-1)', ...
				      'P_PN(100) (mmol N m-2 d-1)', 'P_PN_Std(100) (mmol N m-2 d-1)', ...
				      'Upwelling-adjusted P_PN(100) (mmol N m-2 d-1)', 'Upwelling-adjusted P_PN_Std(100) (mmol N m-2 d-1)', ...
				      'P_PN(E_z) (mmol N m-2 d-1)', 'P_PN_Std(E_z) (mmol N m-2 d-1)', ...
			   	      'Upwelling-adjusted P_PN(E_z) (mmol N m-2 d-1)', 'Upwelling-adjusted P_PN_Std(E_z) (mmol N m-2 d-1)'}; 

%   get corrected values 
%%% th234 at ppz ::
th234PpzFlux = depthData.('P_Th(E_z) (dpm m-2 d-1)'); 
uncertTh234PpzFlux = depthData.('P_Th_Std(E_z) (dpm m-2 d-1)'); 
th234PpzFlux(~isnan(depthData.('Upwelling-adjusted P_Th(E_z) (dpm m-2 d-1)'))) = depthData.('Upwelling-adjusted P_Th(E_z) (dpm m-2 d-1)')(~isnan(depthData.('Upwelling-adjusted P_Th(E_z) (dpm m-2 d-1)'))); 
uncertTh234PpzFlux(~isnan(depthData.('Upwelling-adjusted P_Th_Std(E_z) (dpm m-2 d-1)'))) = depthData.('Upwelling-adjusted P_Th_Std(E_z) (dpm m-2 d-1)')(~isnan(depthData.('Upwelling-adjusted P_Th_Std(E_z) (dpm m-2 d-1)'))); 

%%% th234 at 100 m below ppz ::
th234100PpzFlux = depthData.('P_Th(E_z+100) (dpm m-2 d-1)'); 
uncertTh234100PpzFlux = depthData.('P_Th_Std(E_z+100) (dpm m-2 d-1)'); 
th234100PpzFlux(~isnan(depthData.('Upwelling-adjusted P_Th(E_z+100) (dpm m-2 d-1)'))) = depthData.('Upwelling-adjusted P_Th(E_z+100) (dpm m-2 d-1)')(~isnan(depthData.('Upwelling-adjusted P_Th(E_z+100) (dpm m-2 d-1)'))); 
uncertTh234100PpzFlux(~isnan(depthData.('Upwelling-adjusted P_Th_Std(E_z+100) (dpm m-2 d-1)'))) = depthData.('Upwelling-adjusted P_Th_Std(E_z+100) (dpm m-2 d-1)')(~isnan(depthData.('Upwelling-adjusted P_Th_Std(E_z+100) (dpm m-2 d-1)'))); 

%%% poc at ppz ::
pocPpzFlux = depthData.('P_POC(E_z) (mmol C m-2 d-1)'); 
uncertPocPpzFlux = depthData.('P_POC_Std(E_z) (mmol C m-2 d-1)'); 
pocPpzFlux(~isnan(depthData.('Upwelling-adjusted P_POC(E_z) (mmol C m-2 d-1)'))) = depthData.('Upwelling-adjusted P_POC(E_z) (mmol C m-2 d-1)')(~isnan(depthData.('Upwelling-adjusted P_POC(E_z) (mmol C m-2 d-1)'))); 
uncertPocPpzFlux(~isnan(depthData.('Upwelling-adjusted P_POC_Std(E_z) (mmol C m-2 d-1)'))) = depthData.('Upwelling-adjusted P_POC_Std(E_z) (mmol C m-2 d-1)')(~isnan(depthData.('Upwelling-adjusted P_POC_Std(E_z) (mmol C m-2 d-1)'))); 

%%% poc at 100 m below ppz ::
poc100PpzFlux = depthData.('P_POC(E_z+100) (mmol C m-2 d-1)'); 
uncertPoc100PpzFlux = depthData.('P_POC_Std(E_z+100) (mmol C m-2 d-1)'); 
poc100PpzFlux(~isnan(depthData.('Upwelling-adjusted P_POC(E_z+100) (mmol C m-2 d-1)'))) = depthData.('Upwelling-adjusted P_POC(E_z+100) (mmol C m-2 d-1)')(~isnan(depthData.('Upwelling-adjusted P_POC(E_z+100) (mmol C m-2 d-1)'))); 
% uncertPoc100PpzFlux(~isnan(depthData.('Upwelling-adjusted P_POC_Std(E_z+100) (mmol C m-2 d-1)'))) = depthData.('Upwelling-adjusted P_POC(E_z+100) (mmol C m-2 d-1)')(~isnan(depthData.('Upwelling-adjusted P_POC(E_z+100) (mmol C m-2 d-1)')));  % bug: RHS uses flux column (P_POC) instead of uncertainty column (P_POC_Std)
uncertPoc100PpzFlux(~isnan(depthData.('Upwelling-adjusted P_POC_Std(E_z+100) (mmol C m-2 d-1)'))) = depthData.('Upwelling-adjusted P_POC_Std(E_z+100) (mmol C m-2 d-1)')(~isnan(depthData.('Upwelling-adjusted P_POC_Std(E_z+100) (mmol C m-2 d-1)'))); 

%%% npp standard error (35-day trailing mean, matched to Th-234 integration window) ::
% uncertNpp = sqrt(mean((gp15_npp.mmolC - mean(gp15_npp.mmolC)) .^ 2) ./ length(gp15_npp.mmolC));
% m1: kriging SD alone — underestimates CBPM algorithmic error in well-sampled regions
% uncertNpp = gp15_npp.stdev;
cbpmFrac  = 0.30;   % 30% algorithmic uncertainty (Behrenfeld et al.; Saba et al. 2011)
uncertNpp = sqrt((gp15_npp.stdev_35d .^ 2) + (cbpmFrac .* gp15_npp.mmolC_35d) .^ 2);

%   calculate metrics
%%% calculate T100 ::
depthData.T100 = poc100PpzFlux ./ pocPpzFlux;
% H2: T100 = F_POC(PPZ+100)/F_POC(PPZ) = F_Th(PPZ+100)*R / (F_Th(PPZ)*R). In the delta
%     method for this ratio, the R terms cancel exactly when both depths share the same
%     deep-zone ratio (exact covariance algebra: the sigma_R^2 contributions sum to zero).
%     The uncertainty therefore reduces to the Th-234 nested-covariance formula:
%       Var[T100] = (Var_Th(PPZ+100) - T100*(2-T100)*Var_Th(PPZ)) / F_Th(PPZ)^2
%     Using POC variances in the same expression (incorrect): inflated by sigma_R^2, giving
%     a negative argument -> clamped to 0 by max -> sigma_T100 = 0 (wrong).
%     Adding 2*T100^2*(sigma_R/R)^2 on top (Option B): overcounts, because that term is zero.
% depthData.uncertT100 = sqrt(max(uncertPoc100PpzFlux .^ 2 - depthData.T100 .* (2 - depthData.T100) .* uncertPocPpzFlux .^ 2, 0)) ./ abs(pocPpzFlux);
depthData.uncertT100 = sqrt(max(uncertTh234100PpzFlux .^ 2 - depthData.T100 .* (2 - depthData.T100) .* uncertTh234PpzFlux .^ 2, 0)) ./ abs(th234PpzFlux);

%%% calculate R100 ::
depthData.R100 = th234100PpzFlux - th234PpzFlux;
% S1: Cov[F(PPZ+100), F(PPZ)] = Var[F(PPZ)] for nested cumulative sums (shared layers 1..PPZ)
% Var[R100] = Var[F(PPZ+100)] + Var[F(PPZ)] - 2*Cov = Var[F(PPZ+100)] - Var[F(PPZ)]
% Previous formula double-counted shared-layer variance and overstated sigma(R100).
% depthData.uncertR100 = sqrt((uncertTh234100PpzFlux .^ 2) + (uncertTh234PpzFlux .^ 2));
depthData.uncertR100 = sqrt(max(uncertTh234100PpzFlux .^ 2 - uncertTh234PpzFlux .^ 2, 0));

%%% calculate Ez Ratio using 35-day trailing mean NPP to match Th-234 integration window ::
depthData.EzRatio = pocPpzFlux ./ gp15_npp.mmolC_35d;
% depthData.uncertEzRatio = abs(depthData.EzRatio) .* sqrt(((uncertPoc100PpzFlux ./ poc100PpzFlux) .^ 2) + ((uncertNpp ./ gp15_npp.mmolC_35d) .^ 2));  % bug #8: EzRatio = pocPpzFlux/NPP so numerator uncertainty must use pocPpzFlux, not poc100PpzFlux (PPZ+100 m)
depthData.uncertEzRatio = abs(depthData.EzRatio) .* sqrt(((uncertPocPpzFlux ./ pocPpzFlux) .^ 2) + ((uncertNpp ./ gp15_npp.mmolC_35d) .^ 2));

%   regional averaging ::
%%% loop through all regions ::
for i = 1 : 1 : numRegions

	%   get data ::
	%%% get stations in region ::
	regStatData = gp15_stations((gp15_stations.stationNo >= regions.first(i)) & (gp15_stations.stationNo <= regions.last(i)), :);

	%%% get mld ::
	regMld = gp15_mld((gp15_mld.('Station No') >= regions.first(i)) & (gp15_mld.('Station No') <= regions.last(i)), :);
	
	%%% get npp ::
	regNpp = gp15_npp((gp15_npp.Station >= regions.first(i)) & (gp15_npp.Station <= regions.last(i)), :);

	%%% get particulate data ::
	regPocLsf = ratioPocLsf(i, :);  	
	regPnLsf = ratioPnLsf(i, :);  	
	regPocSsf = ratioPocSsf(i, :);  	
	regPnSsf = ratioPnSsf(i, :);  	

	%%% flux data and metrics ::
	regFlux = depthData((depthData.Station >= regions.first(i)) & (depthData.Station <= regions.last(i)), :);	

	%   perform mean and std calculations ::
	latRange = round([max(regStatData.latitude), min(regStatData.latitude)], 0); 
	stationRange = [min(regStatData.stationNo), max(regStatData.stationNo)]; 
	% S2: nanmean/nanstd deprecated; replaced with mean/std(..., 'omitnan')
	meanMld = round([mean(regMld.('MLD JAK'), 'omitnan'), std(regMld.('MLD JAK'), 0, 'omitnan')], 0);
	meanPpz = round([mean(regStatData.depthPpz, 'omitnan'), std(regStatData.depthPpz, 0, 'omitnan')], 0);
	meanNpp = round([mean(regNpp.mmolC_35d, 'omitnan'), std(regNpp.mmolC_35d, 0, 'omitnan')], 0);
	meanSsLsfCtoTh = round([regPocLsf.ratioSubsurface, regPocLsf.ratioSubsurfaceStd], 1); 
	meanSsSsfCtoTh = round([regPocSsf.ratioSubsurface, regPocSsf.ratioSubsurfaceStd], 1);
	meanSsLsfNtoTh = round([regPnLsf.ratioSubsurface, regPnLsf.ratioSubsurfaceStd], 2);
	meanSsSsfNtoTh = round([regPnSsf.ratioSubsurface, regPnSsf.ratioSubsurfaceStd], 2);
	% S2: nanmean/nanstd deprecated; replaced with mean/std(..., 'omitnan')
	meanPThEz = round([mean(regFlux.('P_Th(E_z) (dpm m-2 d-1)'), 'omitnan'), std(regFlux.('P_Th(E_z) (dpm m-2 d-1)'), 0, 'omitnan')], -2);
	meanPTh100 = round([mean(regFlux.('P_Th(100) (dpm m-2 d-1)'), 'omitnan'), std(regFlux.('P_Th(100) (dpm m-2 d-1)'), 0, 'omitnan')], -2);
	meanPPocEz = round([mean(regFlux.('P_POC(E_z) (mmol C m-2 d-1)'), 'omitnan'), std(regFlux.('P_POC(E_z) (mmol C m-2 d-1)'), 0, 'omitnan')], 2);
	meanPPoc100 = round([mean(regFlux.('P_POC(100) (mmol C m-2 d-1)'), 'omitnan'), std(regFlux.('P_POC(100) (mmol C m-2 d-1)'), 0, 'omitnan')], 2);
	meanPPnEz = round([mean(regFlux.('P_PN(E_z) (mmol N m-2 d-1)'), 'omitnan'), std(regFlux.('P_PN(E_z) (mmol N m-2 d-1)'), 0, 'omitnan')], 2);
	meanPPn100 = round([mean(regFlux.('P_PN(100) (mmol N m-2 d-1)'), 'omitnan'), std(regFlux.('P_PN(100) (mmol N m-2 d-1)'), 0, 'omitnan')], 2);
	meanT100 = round([mean(regFlux.T100, 'omitnan'), std(regFlux.T100, 0, 'omitnan')], 1);
	meanR100 = round([mean(regFlux.R100, 'omitnan'), std(regFlux.R100, 0, 'omitnan')], 1);
	meanEzRatio = round([mean(regFlux.EzRatio, 'omitnan'), std(regFlux.EzRatio, 0, 'omitnan')], 2);

	%   store data ::
	regionData{i} = {latRange, stationRange, meanMld, meanPpz, meanNpp, ...
	                 meanSsLsfCtoTh, meanSsSsfCtoTh, meanSsLsfNtoTh, meanSsSsfNtoTh, ...
			 meanPThEz, meanPTh100, meanPPocEz, meanPPoc100, meanPPnEz, meanPPn100, ...
			 meanT100, meanR100, meanEzRatio};

end

%   make tables ::
%%% regional data ::
for i = 1 : 1 : numRegions
	dat = cell2table(regionData{i}'); 
	dat.Properties.VariableNames = {regionNames{i}}; 
	dat.Properties.RowNames = {'Latitudinal Range', 'Stations', 'MLD', 'PPZ', 'NPP', ...
	                           'LSF C:Th', 'SSF C:Th', 'LSF N:Th', 'SSF N:Th', ...
				   'P_Th(E_z)', 'P_Th(100)', 'P_POC(E_z)', 'P_POC(100)', ...
				   'P_PN(E_z)', 'P_PN(100)', 'T_100', 'R_100', 'Ez-Ratio'};
	if i == 1
		regionalData = dat; 
	else
		regionalData = [regionalData, dat]; 
	end
end

%%% table 1 ::
table1 = depthData(:, 1:12); 

%%% table 2 ::
table2 = regionalData([1:16, 18], :); 

%   write data ::
%%% .mat ::
save([sim_output_basepath 'gp15Model/modelOutput/depthData.mat'], 'depthData'); 
save([sim_output_basepath 'gp15Model/modelOutput/regionData.mat'], 'regionData'); 
save([sim_output_basepath 'gp15Model/tables/table1.mat'], 'table1'); 
save([sim_output_basepath 'gp15Model/tables/table2.mat'], 'table2'); 

%%% .xlsx ::
writetable(depthData, [sim_output_basepath 'gp15Model/modelOutput/depthData.xlsx']); 
writetable(regionalData, [sim_output_basepath 'gp15Model/modelOutput/regionData.xlsx'], 'writeRowNames', 1); 
writetable(table1, [sim_output_basepath 'gp15Model/tables/table1.xlsx']); 
writetable(table2, [sim_output_basepath 'gp15Model/tables/table2.xlsx'], 'writeRowNames', 1); 

%%  calculate text statistics
%   load observational data ::
load([pro_output_basepath 'doQc/gp15/gp15_obs.mat'], 'gp15_obs');
load([pro_output_basepath 'doQc/gp15/gp15_obsNoQc.mat'], 'gp15_obsNoQc');

%   3.2 ::
%%% range ::
textStats.results32.th234Range = round([min(gp15_obs.th234), max(gp15_obs.th234)], 1); 

%%% ppz range ::
textStats.results32.ppzRange = round([min(gp15_stations.depthPpz), max(gp15_stations.depthPpz)], 0);
textStats.results32.ppzRange = [gp15_stations.stationNo(round(gp15_stations.depthPpz, 0) == textStats.results32.ppzRange(1)), gp15_stations.stationNo(round(gp15_stations.depthPpz, 0) == textStats.results32.ppzRange(2)); textStats.results32.ppzRange];

%%% lsf ::
lsfFrac = gp15_obs.th234PocLarge ./ gp15_obs.th234;
lsfFrac(lsfFrac > 1) = NaN; 
idx100 = (gp15_obs.depth >= 0) & (gp15_obs.depth < 100); 
idx200 = (gp15_obs.depth >= 100) & (gp15_obs.depth < 200); 
% S2: nanmean deprecated
textStats.results32.lsf234Th = round([1, min(lsfFrac(idx100)), max(lsfFrac(idx100)), mean(lsfFrac(idx100), 'omitnan'); ...
				      2, min(lsfFrac(idx200)), max(lsfFrac(idx200)), mean(lsfFrac(idx200), 'omitnan')] * 100, 0);

%%% ssf ::
ssfFrac = gp15_obs.th234PocSmall ./ gp15_obs.th234;
ssfFrac(ssfFrac > 1) = NaN; 
% S2: nanmean deprecated
textStats.results32.ssf234Th = round([1, min(ssfFrac(idx100)), max(ssfFrac(idx100)), mean(ssfFrac(idx100), 'omitnan'); ...
				      2, min(ssfFrac(idx200)), max(ssfFrac(idx200)), mean(ssfFrac(idx200), 'omitnan')] * 100, 0);

%   3.3 ::
%%% npp (35-day trailing mean, matched to Th-234 integration window) ::
textStats.results33.nppRange = round([min(gp15_npp.mmolC_35d), max(gp15_npp.mmolC_35d)], 0);
textStats.results33.nppStatRange.minNpp = [gp15_npp.Station(round(gp15_npp.mmolC_35d, 0) == textStats.results33.nppRange(1)), ...
				           ones(size(gp15_npp.Station(round(gp15_npp.mmolC_35d, 0) == textStats.results33.nppRange(1)))) * textStats.results33.nppRange(1)];
textStats.results33.nppStatRange.maxNpp = [gp15_npp.Station(round(gp15_npp.mmolC_35d, 0) == textStats.results33.nppRange(2)), ...
				           ones(size(gp15_npp.Station(round(gp15_npp.mmolC_35d, 0) == textStats.results33.nppRange(2)))) * textStats.results33.nppRange(2)];

%%% lsf poc:234th ratio anomalous ::
lsfPocRatioSort = sort(gp15_obsNoQc.pocTh234RatioLarge(~isnan(gp15_obsNoQc.pocTh234RatioLarge))); 
textStats.results33.lsfPocRatioAnamalous = round([lsfPocRatioSort(end-1), lsfPocRatioSort(end)], 0); 
textStats.results33.lsfPocRatioAnamalous = [gp15_obsNoQc.stationNo(round(gp15_obsNoQc.pocTh234RatioLarge, 0) == textStats.results33.lsfPocRatioAnamalous(1)), gp15_obsNoQc.stationNo(round(gp15_obsNoQc.pocTh234RatioLarge, 0) == textStats.results33.lsfPocRatioAnamalous(2)); textStats.results33.lsfPocRatioAnamalous];

%%% lsf poc:234th ratio ::
lsfPocRatio = sort(gp15_obs.pocTh234RatioLarge(~isnan(gp15_obs.pocTh234RatioLarge)));
maxPocRatio = round(max(lsfPocRatio), 0); 
textStats.results33.lsfPocRatio = [gp15_obs.stationNo(round(gp15_obs.pocTh234RatioLarge, 0) == maxPocRatio), maxPocRatio];

%%% negative flux ::
idxNegFlux = (gp15_fluxPpz.th234FluxCumul1d < 0);
textStats.results33.negativeEzFlux = [gp15_fluxPpz.stationNo(idxNegFlux), round(gp15_fluxPpz.th234FluxCumul1d(idxNegFlux), -1), round(gp15_fluxPpz.uncert1d(idxNegFlux), -1)];

%%% correction by regression ::
pocFluxCompare = [depthData.('P_Th(E_z) (dpm m-2 d-1)') .* gp15_fluxPpz.pocTh234RatioLarge, depthData.('P_Th(E_z) (dpm m-2 d-1)') .* gp15_fluxPpz.pocTh234RatioRegress];
textStats.results33.pocRegressCorrect = round(nanmedian(abs(pocFluxCompare(:, 1) - pocFluxCompare(:, 2)) ./ abs(pocFluxCompare(:, 2))), 2); 

%%% lsf pn:234th ratio ::
lsfPnRatio = sort(gp15_obs.pnTh234RatioLarge(~isnan(gp15_obs.pnTh234RatioLarge)));
maxPnRatio = round(lsfPnRatio(end), 0); 
textStats.results33.lsfPnRatioAnamalous = [gp15_obs.stationNo(round(gp15_obs.pnTh234RatioLarge, 0) == maxPnRatio), maxPnRatio];
maxPnRatio = round(lsfPnRatio(end-1), 1); 
textStats.results33.lsfPnRatio = [gp15_obs.stationNo(round(gp15_obs.pnTh234RatioLarge, 1) == maxPnRatio), maxPnRatio];

%%% ssf assumption ::
ssfPocFluxCompare100 = [depthData.('P_Th(100) (dpm m-2 d-1)') .* gp15_flux100.pocTh234RatioLarge, depthData.('P_Th(100) (dpm m-2 d-1)') .* gp15_flux100.pocTh234RatioSmall];
ssfPocFluxCompareEz = [depthData.('P_Th(E_z) (dpm m-2 d-1)') .* gp15_fluxPpz.pocTh234RatioLarge, depthData.('P_Th(E_z) (dpm m-2 d-1)') .* gp15_fluxPpz.pocTh234RatioSmall];

%   3.4 ::
%%% highest flux ::
maxFlux = unique(round(sort(gp15_fluxPpz.th234FluxCumul1d(gp15_fluxPpz.th234FluxCumul1d > 0)), -2)); 
maxFlux = maxFlux(end-1:end); 
textStats.results34.maxFlux = [gp15_fluxPpz.stationNo(round(gp15_fluxPpz.th234FluxCumul1d, -2) == maxFlux(1)), gp15_fluxPpz.stationNo(round(gp15_fluxPpz.th234FluxCumul1d, -2) == maxFlux(2)); maxFlux'];

%%% lowest flux ::
minFlux = round(sort(gp15_fluxPpz.th234FluxCumul1d(gp15_fluxPpz.th234FluxCumul1d > 0)), -1); 
minFlux = minFlux(1:2); 
textStats.results34.minFlux = [gp15_fluxPpz.stationNo(round(gp15_fluxPpz.th234FluxCumul1d, -1) == minFlux(1)), gp15_fluxPpz.stationNo(round(gp15_fluxPpz.th234FluxCumul1d, -1) == minFlux(2)); minFlux'];

%   4.1 ::
%%% upwelling correction at 100 m ::
% upwellCorrectStats = gp15_flux100.stationNo;  % dead variable: assigned but never read
upwellCorrect = (gp15_flux100.(['th234FluxCumul2d_' fluxProduct]) - gp15_flux100.th234FluxCumul1d) ./ gp15_flux100.th234FluxCumul1d;
% S2: nanmean/nanstd deprecated
textStats.results41.upwellCorrectRange100 = round([mean(upwellCorrect, 'omitnan'), std(upwellCorrect, 0, 'omitnan'), min(upwellCorrect), max(upwellCorrect)] * 100, 0);

%%% upwelling correction at ppz ::
% upwellCorrectStats = gp15_fluxPpz.stationNo;  % dead variable: assigned but never read
upwellCorrect = (gp15_fluxPpz.(['th234FluxCumul2d_' fluxProduct]) - gp15_fluxPpz.th234FluxCumul1d) ./ gp15_fluxPpz.th234FluxCumul1d;
% S2: nanmean/nanstd deprecated
textStats.results41.upwellCorrectRangePpz = round([mean(upwellCorrect, 'omitnan'), std(upwellCorrect, 0, 'omitnan'), min(upwellCorrect), max(upwellCorrect)] * 100, 0);

%   4.4 ::
%%% ez ratios ::
textStats.results44.ezRatios = round(max(depthData.EzRatio), 2); 

%%% npp at 100 below ez ::
npp100Ppz = depthData.EzRatio .* depthData.T100; 
textStats.results44.npp100Ppz = max(npp100Ppz); 

%%% npp eqpac ::
nppEqPac = gp15_npp.mmolC_35d((gp15_npp.Station >= regions.first(3)) & (gp15_npp.Station <= regions.last(3)));
textStats.results44.maxNppEqPac = round(max(nppEqPac), 0); 

%%  save data
save([sim_output_basepath 'gp15Model/modelOutput/textStats.mat'], 'textStats');

%%  end subroutine
