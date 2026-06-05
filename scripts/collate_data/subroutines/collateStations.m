%%  collect data per stations
%   n.b.: just take gp15_obs and add in upwelling and vertical gradient, using depth as a key per station.
%         also, cut off at bottom depth. make sure to keep unit conversion, though.
%   set coefficients ::
L2M = 1000;
UMOL2MMOL = 1000;
SEC2DAY = 60 * 60 * 24;

%   set max flux depth ::
BTMDEPTH = 400;  % m, THIS IS VERY IMPORTANT

%   loop through all stations ::
for iStat = 1 : 1 : NUMSTAT

	%   get station number ::
	statNo = gp15_stations.stationNo(iStat);

	%   get station data ::
	statData = gp15_obs((gp15_obs.stationNo == statNo) & (gp15_obs.depth <= BTMDEPTH), :);

	%   get upwelling data ::
	idxStat = (gp15_w.ecco.stationNo == statNo) & (gp15_w.ecco.depth <= BTMDEPTH);
	statUpwell.depth = table2array(gp15_w.ecco(idxStat, 'depth'));
	statUpwell.w_ecco = table2array(gp15_w.ecco(idxStat, 'w'));
	statUpwell.wErr_ecco = table2array(gp15_w.ecco(idxStat, 'wErr'));
	for dataProduct = {'cglo', 'foam', 'glor', 'oras', 'grep'}
		% E4: eval removed; use dynamic struct field access
		% modelName = ['gp15_w.' dataProduct{1}];
		% eval(['statUpwell.w_'    dataProduct{1} ' = table2array(' modelName '(idxStat, ''w''));']);
		% eval(['statUpwell.wErr_' dataProduct{1} ' = table2array(' modelName '(idxStat, ''wErr''));']);
		statUpwell.(['w_' dataProduct{1}])    = table2array(gp15_w.(dataProduct{1})(idxStat, 'w'));
		statUpwell.(['wErr_' dataProduct{1}]) = table2array(gp15_w.(dataProduct{1})(idxStat, 'wErr'));
	end
	statUpwell = struct2table(statUpwell);

	%   get vertical gradient ::
	statVertGrad = gp15_grad((gp15_grad.stationNo == statNo) & (gp15_grad.depth <= BTMDEPTH), {'depth', 'vertGrad', 'vertGradError'});

	%   get regression ratio ::
	statRatioRegressPoc = regressPocRatio((regressPocRatio.stationNo == statNo) & (regressPocRatio.depth <= BTMDEPTH), {'depth', 'pocTh234RatioRegress', 'uncertPocTh234RatioRegress'});
	statRatioRegressPn = regressPnRatio((regressPnRatio.stationNo == statNo) & (regressPnRatio.depth <= BTMDEPTH), {'depth', 'pnTh234RatioRegress', 'uncertPnTh234RatioRegress'});

	%   get regions ::
	statRegions = gp15Regions((gp15Regions.stationNo == statNo) & (gp15Regions.depth <= BTMDEPTH), {'depth', 'region'});

	%   assert unique depths — join() silently multiplies rows if any table has duplicates ::
	assert(numel(unique(statData.depth))            == height(statData),            'Station %d: duplicate depths in gp15_obs.',        statNo);
	assert(numel(unique(statUpwell.depth))          == height(statUpwell),          'Station %d: duplicate depths in gp15_w.',          statNo);
	assert(numel(unique(statVertGrad.depth))        == height(statVertGrad),        'Station %d: duplicate depths in gp15_grad.',       statNo);
	assert(numel(unique(statRatioRegressPoc.depth)) == height(statRatioRegressPoc), 'Station %d: duplicate depths in regressPocRatio.', statNo);
	assert(numel(unique(statRatioRegressPn.depth))  == height(statRatioRegressPn),  'Station %d: duplicate depths in regressPnRatio.',  statNo);
	assert(numel(unique(statRegions.depth))         == height(statRegions),         'Station %d: duplicate depths in gp15Regions.',     statNo);

	%   correct units ::
	%%% radionuclide ::
	% E3: column index {9:12} silently corrupts if gp15_obs gains a column; use named access
	% statData{:, 9:12} = statData{:, 9:12} * L2M;
	statData.th234       = statData.th234       * L2M;
	statData.uncertTh234 = statData.uncertTh234 * L2M;
	statData.u238        = statData.u238        * L2M;
	statData.uncertU238  = statData.uncertU238  * L2M;

	%%% particulate ::
	% E3: column index {13:end} silently corrupts if gp15_obs gains a column; use named access
	% statData{:, 13:end} = statData{:, 13:end} / UMOL2MMOL;
	statData.pocLarge                 = statData.pocLarge                 / UMOL2MMOL;
	statData.uncertPocLarge           = statData.uncertPocLarge           / UMOL2MMOL;
	statData.th234PocLarge            = statData.th234PocLarge            / UMOL2MMOL;
	statData.uncertTh234PocLarge      = statData.uncertTh234PocLarge      / UMOL2MMOL;
	statData.pocTh234RatioLarge       = statData.pocTh234RatioLarge       / UMOL2MMOL;
	statData.uncertPocTh234RatioLarge = statData.uncertPocTh234RatioLarge / UMOL2MMOL;
	statData.pocSmall                 = statData.pocSmall                 / UMOL2MMOL;
	statData.uncertPocSmall           = statData.uncertPocSmall           / UMOL2MMOL;
	statData.th234PocSmall            = statData.th234PocSmall            / UMOL2MMOL;
	statData.uncertTh234PocSmall      = statData.uncertTh234PocSmall      / UMOL2MMOL;
	statData.pocTh234RatioSmall       = statData.pocTh234RatioSmall       / UMOL2MMOL;
	statData.uncertPocTh234RatioSmall = statData.uncertPocTh234RatioSmall / UMOL2MMOL;
	statData.pnLarge                  = statData.pnLarge                  / UMOL2MMOL;
	statData.uncertPnLarge            = statData.uncertPnLarge            / UMOL2MMOL;
	statData.pnTh234RatioLarge        = statData.pnTh234RatioLarge        / UMOL2MMOL;
	statData.uncertPnTh234RatioLarge  = statData.uncertPnTh234RatioLarge  / UMOL2MMOL;
	statData.pnSmall                  = statData.pnSmall                  / UMOL2MMOL;
	statData.uncertPnSmall            = statData.uncertPnSmall            / UMOL2MMOL;
	statData.pnTh234RatioSmall        = statData.pnTh234RatioSmall        / UMOL2MMOL;
	statData.uncertPnTh234RatioSmall  = statData.uncertPnTh234RatioSmall  / UMOL2MMOL;

	%%% velocity ::
	% E3: column index {2:end} silently corrupts if statUpwell gains a column; use named access via loop
	% statUpwell{:, 2:end} = statUpwell{:, 2:end} * SEC2DAY;
	for dataProduct = {'ecco', 'cglo', 'foam', 'glor', 'oras', 'grep'}
		dp = dataProduct{1};
		statUpwell.(['w_' dp])    = statUpwell.(['w_' dp])    * SEC2DAY;
		statUpwell.(['wErr_' dp]) = statUpwell.(['wErr_' dp]) * SEC2DAY;
	end

	%%% ensemble mean velocity (added after per-product unit conversion) ::
	%   σ_mean = sqrt(Σ σ_k²) / n  — correct propagation for uncertainty of a mean of n independent quantities
	statUpwell.w_meanFull      = mean([statUpwell.w_ecco, statUpwell.w_cglo, statUpwell.w_foam, statUpwell.w_glor, statUpwell.w_oras], 2, 'omitnan');
	statUpwell.wErr_meanFull   = sqrt(statUpwell.wErr_ecco.^2 + statUpwell.wErr_cglo.^2 + statUpwell.wErr_foam.^2 + statUpwell.wErr_glor.^2 + statUpwell.wErr_oras.^2) / 5;
	statUpwell.w_meanNoFOAM    = mean([statUpwell.w_ecco, statUpwell.w_cglo, statUpwell.w_glor, statUpwell.w_oras], 2, 'omitnan');
	statUpwell.wErr_meanNoFOAM = sqrt(statUpwell.wErr_ecco.^2 + statUpwell.wErr_cglo.^2 + statUpwell.wErr_glor.^2 + statUpwell.wErr_oras.^2) / 4;

	%%% gradient values ::
	% E3: column index {2:3} silently corrupts; use named access
	% statVertGrad{:, 2:3} = statVertGrad{:, 2:3} * L2M;
	statVertGrad.vertGrad      = statVertGrad.vertGrad      * L2M;
	statVertGrad.vertGradError = statVertGrad.vertGradError * L2M;

	%%% regression ratio ::
	% E3: column index {2:3} silently corrupts; use named access
	% statRatioRegressPoc{:, 2:3} = statRatioRegressPoc{:, 2:3} / UMOL2MMOL;
	% statRatioRegressPn{:, 2:3} = statRatioRegressPn{:, 2:3} / UMOL2MMOL;
	statRatioRegressPoc.pocTh234RatioRegress       = statRatioRegressPoc.pocTh234RatioRegress       / UMOL2MMOL;
	statRatioRegressPoc.uncertPocTh234RatioRegress = statRatioRegressPoc.uncertPocTh234RatioRegress / UMOL2MMOL;
	statRatioRegressPn.pnTh234RatioRegress         = statRatioRegressPn.pnTh234RatioRegress         / UMOL2MMOL;
	statRatioRegressPn.uncertPnTh234RatioRegress   = statRatioRegressPn.uncertPnTh234RatioRegress   / UMOL2MMOL;
	statRatioRegress = join(statRatioRegressPoc, statRatioRegressPn, 'keys', 'depth');

	%   make velocity zero where gradient zero ::
	% E3: column index {2} and {2:end} silently corrupts; use named access
	% statUpwell{statVertGrad{:, 2} == 0, 2:end} = 0;
	% S2: assert depth grids match before positional index is applied ::
	assert(isequal(statUpwell.depth, statVertGrad.depth), 'Station %d: upwelling/gradient depth grids differ.', statNo);
	zeroGradIdx = statVertGrad.vertGrad == 0;
	for dataProduct = {'ecco', 'cglo', 'foam', 'glor', 'oras', 'grep', 'meanFull', 'meanNoFOAM'}
		dp = dataProduct{1};
		statUpwell.(['w_' dp])(zeroGradIdx)    = 0;
		statUpwell.(['wErr_' dp])(zeroGradIdx) = 0;
	end

	%   make array ::
	gp15_inputsAdd = join(join(join(join(statRegions, statData, 'keys', 'depth'), statUpwell, 'keys', 'depth'), statVertGrad, 'keys', 'depth'), statRatioRegress, 'keys', 'depth');

	if iStat == 1

		gp15_inputs = gp15_inputsAdd;

	else

		gp15_inputs = vertcat(gp15_inputs, gp15_inputsAdd);

	end

	%   clear ::
	clear('statUpwell');

end

%%  end subroutine
