% see plotQcPoc.m for diagnostic figures

%   set flier values ::
FLIERDATA = 30; 

%   initialize ::
%%% make poc lsf ::
len_data = size(gp15_obs, 1); 
gp15_obs.pocLarge = NaN(len_data, 1);
gp15_obs.uncertPocLarge = NaN(len_data, 1);
gp15_obs.th234PocLarge = NaN(len_data, 1);
gp15_obs.uncertTh234PocLarge = NaN(len_data, 1);
gp15_obs.pocTh234RatioLarge = NaN(len_data, 1);
gp15_obs.uncertPocTh234RatioLarge = NaN(len_data, 1);

%%% make poc ssf ::
gp15_obs.pocSmall = NaN(len_data, 1);
gp15_obs.uncertPocSmall = NaN(len_data, 1);
gp15_obs.th234PocSmall = NaN(len_data, 1);
gp15_obs.uncertTh234PocSmall = NaN(len_data, 1);
gp15_obs.pocTh234RatioSmall = NaN(len_data, 1);
gp15_obs.uncertPocTh234RatioSmall = NaN(len_data, 1);

%   loop through all stations ::
for iStat = 1 : 1 : NUMSTAT

	%   get station number ::
	statNo = gp15_stations.stationNo(iStat);

	%   get station data ::
	statData = gp15_obs(gp15_obs.stationNo == statNo, :);
	statDataRadioPoc = gp15_part(gp15_part.stationNo == statNo, :);
	statDataPoc = gp15_particles(gp15_particles.station_num == statNo, :);  

	%    get depths ::
	xq = statData.depth;

	%    quality control negative values ::
	statDataRadioPoc.th234PartLsf(statDataRadioPoc.th234PartLsf < 0) = NaN;
	statDataRadioPoc.th234PartSsf(statDataRadioPoc.th234PartSsf < 0) = NaN;

	%    get sample data ::
	XRadioPocLarge = table2array(sortrows(rmmissing(statDataRadioPoc(:, {'depth', 'th234PartLsf', 'uncertTh234PartLsf'})), 'depth'));
	XRadioPocSmall = table2array(sortrows(rmmissing(statDataRadioPoc(:, {'depth', 'th234PartSsf', 'uncertTh234PartSsf'})), 'depth'));
	XPocLarge = table2array(sortrows(rmmissing(statDataPoc(:, {'depth_m', 'POC_LPT_uM', 'POC_LPT_uM_std'})), 'depth_m'));
	XPocSmall = table2array(sortrows(rmmissing(statDataPoc(:, {'depth_m', 'POC_SPT_uM', 'POC_SPT_uM_std'})), 'depth_m'));

	%    deal with empty arrays ::
	if isempty(XRadioPocLarge)
		XRadioPocLarge = NaN(size(statData, 1), 3);
	end
	if isempty(XRadioPocSmall)
		XRadioPocSmall = NaN(size(statData, 1), 3);
	end
	if isempty(XPocLarge)
		XPocLarge = NaN(size(statData, 1), 3);
	end
	if isempty(XPocSmall)
		XPocSmall = NaN(size(statData, 1), 3);
	end

	%    interpolate ::
	[th234PocLarge, th234PocLargeStd] = interp1dError(XRadioPocLarge(:, 1), XRadioPocLarge(:, 2), XRadioPocLarge(:, 3), xq); 
	[th234PocSmall, th234PocSmallStd] = interp1dError(XRadioPocSmall(:, 1), XRadioPocSmall(:, 2), XRadioPocSmall(:, 3), xq); 
	[pocLarge, pocLargeStd] = interp1dError(XPocLarge(:, 1), XPocLarge(:, 2), XPocLarge(:, 3), xq); 
	[pocSmall, pocSmallStd] = interp1dError(XPocSmall(:, 1), XPocSmall(:, 2), XPocSmall(:, 3), xq); 
	
	%    make ratio data :: 
	pocTh234RatioLarge = pocLarge ./ th234PocLarge;
	pocTh234RatioLargeStd = abs(pocTh234RatioLarge) .* sqrt(((pocLargeStd ./ pocLarge) .^ 2) + ((th234PocLargeStd ./ th234PocLarge) .^ 2));
	pocTh234RatioSmall = pocSmall ./ th234PocSmall;
	pocTh234RatioSmallStd = abs(pocTh234RatioSmall) .* sqrt(((pocSmallStd ./ pocSmall) .^ 2) + ((th234PocSmallStd ./ th234PocSmall) .^ 2));
	
	%   put into data array ::
	statData.pocLarge = pocLarge;
	statData.uncertPocLarge = pocLargeStd;
	statData.th234PocLarge = th234PocLarge;
	statData.uncertTh234PocLarge = th234PocLargeStd;
	statData.pocTh234RatioLarge = pocTh234RatioLarge;
	statData.uncertPocTh234RatioLarge = pocTh234RatioLargeStd;
	statData.pocSmall = pocSmall;
	statData.uncertPocSmall = pocSmallStd;
	statData.th234PocSmall = th234PocSmall;
	statData.uncertTh234PocSmall = th234PocSmallStd;
	statData.pocTh234RatioSmall = pocTh234RatioSmall;
	statData.uncertPocTh234RatioSmall = pocTh234RatioSmallStd;
	
	%   put into storage array ::
	if iStat == 1
		gp15_obs_final = statData;
	else
		gp15_obs_final = [gp15_obs_final; statData];
	end

end

%%  hand off data
gp15_obs = gp15_obs_final;
gp15_obsNoQc = gp15_obs_final;

%%  correct flier data
%   remove ::
% gp15_obs.pocTh234RatioLarge(gp15_obs.pocTh234RatioLarge > FLIERDATA) = NaN;   % bug #6: ratio is NaN-ed first, so the uncertainty lines below always evaluate NaN>30 as false and never trigger
% gp15_obs.pocTh234RatioSmall(gp15_obs.pocTh234RatioSmall > FLIERDATA) = NaN;
% gp15_obs.uncertPocTh234RatioLarge(gp15_obs.pocTh234RatioLarge > FLIERDATA) = NaN;
% gp15_obs.uncertPocTh234RatioSmall(gp15_obs.pocTh234RatioSmall > FLIERDATA) = NaN;
flierIdxLarge = gp15_obs.pocTh234RatioLarge > FLIERDATA;
flierIdxSmall = gp15_obs.pocTh234RatioSmall > FLIERDATA;
gp15_obs.pocTh234RatioLarge(flierIdxLarge) = NaN;
gp15_obs.uncertPocTh234RatioLarge(flierIdxLarge) = NaN;
gp15_obs.pocTh234RatioSmall(flierIdxSmall) = NaN;
gp15_obs.uncertPocTh234RatioSmall(flierIdxSmall) = NaN;

%%  clean up
clear('gp15_obs_final');
close('all'); 

%%  end subsubroutine
