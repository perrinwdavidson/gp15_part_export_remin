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
	for dataProduct = {"cglo", "foam", "glor", "oras"}  % grep is the average of the four model ensemble
		modelName = strcat("gp15_w.", dataProduct{1});
		eval(strcat("statUpwell.w_", dataProduct{1}, " = table2array(", modelName, "(idxStat, 'w'));")); 
		eval(strcat("statUpwell.wErr_", dataProduct{1}, " = table2array(", modelName, "(idxStat, 'wErr'));")); 
	end
	statUpwell = struct2table(statUpwell); 

	%   get vertical gradient ::
	statVertGrad = gp15_grad((gp15_grad.stationNo == statNo) & (gp15_grad.depth <= BTMDEPTH), {'depth', 'vertGrad', 'vertGradError'}); 

	%   get regression ratio ::
	statRatioRegressPoc = regressPocRatio((regressPocRatio.stationNo == statNo) & (regressPocRatio.depth <= BTMDEPTH), {'depth', 'pocTh234RatioRegress', 'uncertPocTh234RatioRegress'}); 
	statRatioRegressPn = regressPnRatio((regressPnRatio.stationNo == statNo) & (regressPnRatio.depth <= BTMDEPTH), {'depth', 'pnTh234RatioRegress', 'uncertPnTh234RatioRegress'}); 

	%   get regions ::
	statRegions = gp15Regions((gp15Regions.stationNo == statNo) & (gp15Regions.depth <= BTMDEPTH), {'depth', 'region'}); 

	%   correct units ::
	%%% radionuclide ::
	statData{:, 9:12} = statData{:, 9:12} * L2M; 
	
	%%% particulate ::
	statData{:, 13:end} = statData{:, 13:end} / UMOL2MMOL; 

	%%% velocity ::
	statUpwell{:, 2:end} = statUpwell{:, 2:end} * SEC2DAY; 

	%%% gradient values ::
	statVertGrad{:, 2:3} = statVertGrad{:, 2:3} * L2M;

	%%% regression ratio ::
	statRatioRegressPoc{:, 2:3} = statRatioRegressPoc{:, 2:3} / UMOL2MMOL;
	statRatioRegressPn{:, 2:3} = statRatioRegressPn{:, 2:3} / UMOL2MMOL;
	statRatioRegress = join(statRatioRegressPoc, statRatioRegressPn, 'keys', 'depth'); 

	%   make velocity zero where gradient zero ::
	statUpwell{statVertGrad{:, 2} == 0, 2:end} = 0;

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
