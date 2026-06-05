%%  calculate data
%   get indices ::
[st, ca, Is, Ie, id] = indexstas(gp15_ctd.stationNo, gp15_ctd.castNo);

%   calculate ppz ::
[maxf, maxd, ppz, maxd2, bkg] = findppz(st, Ie, Is, gp15_ctd.ctd_depth, gp15_ctd.ctd_fluorescence);

%   make array ::
ppzCalc = [st, ca, maxd, ppz];

%  quality control the data ::
ppzCalc(ppzCalc(:, 4) == -9999, 4) = NaN; 
ppzCalc((ppzCalc(:, 2) == 0) | (ppzCalc(:, 3) == 0) | (ppzCalc(:, 4) == 0)) = NaN; 
ppzCalc = rmmissing(ppzCalc); 

%  make table ::
ppzCalc = array2table(ppzCalc); 
ppzCalc.Properties.VariableNames = {'stationNo', 'castNo', 'maxFDepth', 'ppzDepth'};

%   only keep certain stations and casts ::
[~, idx_gp15, ~] = intersect(ppzCalc(:, 1:2), gp15_stat_cast); 
gp15_ppz = ppzCalc(idx_gp15, :); 

%   add in coordinates ::
%%% loop through all observations ::
for i = 1 : 1 : size(gp15_ppz, 1)

	% make and store data ::
	gp15_ppz.latitude(i) = mean(gp15_obs.latitude(gp15_obs.stationNo == gp15_ppz.stationNo(i))); 
	gp15_ppz.longitude(i) = mean(gp15_obs.longitude(gp15_obs.stationNo == gp15_ppz.stationNo(i)));
	gp15_ppz.samplingDate(i) = mean(datenum(gp15_obs.date(gp15_obs.stationNo == gp15_ppz.stationNo(i))));

end

%   add in normalization ::
gp15_ppz.normalizationTo100Meters = 100 * ones(size(gp15_ppz.ppzDepth)); 
gp15_ppz.maxDepth400Meters = 400 * ones(size(gp15_ppz.ppzDepth));

%%  end subroutine
