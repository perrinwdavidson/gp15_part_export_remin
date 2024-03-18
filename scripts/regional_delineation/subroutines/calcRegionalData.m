%%  make variables of interest ::
%   preallocate ::
delineate_data_raw = struct();

%   load temperature ::
delineate_data_raw.temperature.stationNo = gp15_ctd.stationNo;
delineate_data_raw.temperature.castNo = gp15_ctd.castNo;
delineate_data_raw.temperature.depth = gp15_ctd.ctd_depth;
delineate_data_raw.temperature.data = gp15_ctd.ctd_temperature;
delineate_data_raw.temperature.flag = gp15_ctd.ctd_temperature_flag;

%   load transmissometery data ::
delineate_data_raw.transmissometer.stationNo = gp15_ctd.stationNo;
delineate_data_raw.transmissometer.castNo = gp15_ctd.castNo;
delineate_data_raw.transmissometer.depth = gp15_ctd.ctd_depth;
delineate_data_raw.transmissometer.data = gp15_ctd.ctd_xmiss;
delineate_data_raw.transmissometer.flag = gp15_ctd.ctd_xmiss_flag;

%   load phosphate data ::
delineate_data_raw.po4.stationNo = gp15_bottles.STNNBR;
delineate_data_raw.po4.castNo = gp15_bottles.CASTNO;
delineate_data_raw.po4.depth = gp15_bottles.DEPTH;
delineate_data_raw.po4.data = gp15_bottles.('PHOSPHATE_D_CONC_BOTTLE::ODF');
delineate_data_raw.po4.flag = gp15_bottles.('PHOSPHATE_D_CONC_BOTTLE::ODF_FLAG_W');

%   load nitrate data ::
delineate_data_raw.no3.stationNo = gp15_bottles.STNNBR;
delineate_data_raw.no3.castNo = gp15_bottles.CASTNO;
delineate_data_raw.no3.depth = gp15_bottles.DEPTH;
delineate_data_raw.no3.data = gp15_bottles.('NITRATE_D_CONC_BOTTLE::ODF');
delineate_data_raw.no3.flag = gp15_bottles.('NITRATE_D_CONC_BOTTLE::ODF_FLAG_W');

%   load pn ::
delineate_data_raw.pn.stationNo = gp15_obs.stationNo;
delineate_data_raw.pn.depth = gp15_obs.depth;
delineate_data_raw.pn.data = gp15_obs.pnLarge + gp15_obs.pnSmall;

%   load hplc ::
delineate_data_raw.hplc.stationNo = gp15_hplc.station;
delineate_data_raw.hplc.data.pico = gp15_hplc.pico;
delineate_data_raw.hplc.data.nano = gp15_hplc.nano;
delineate_data_raw.hplc.data.micro = gp15_hplc.micro;

%%  quality control
%   preallocate ::
temperature = NaN(NUMSTAT, 1);
transmissometer = NaN(NUMSTAT, 1);
phosphate = NaN(NUMSTAT, 1);
nitrate = NaN(NUMSTAT, 1);
particulate_n = NaN(NUMSTAT, 1);

%   temperature ::
for istat = 1 : 1 : NUMSTAT

	% get station and cast numbers ::
	stat = gp15_stations.stationNo(istat);  
	cast = gp15_stations.castNo(istat);  

	% get ppz (to match pigment calculation) ::
	ppz = gp15_stations.depthPpz(istat); 

	% get station data ::
	temp = delineate_data_raw.temperature.data((delineate_data_raw.temperature.stationNo == stat) & ...
						   (delineate_data_raw.temperature.castNo == cast) & ...
						   ((delineate_data_raw.temperature.flag == 1) | (delineate_data_raw.temperature.flag == 2)) & ...
					           (delineate_data_raw.temperature.depth <= ppz) & ...
						   ~isnan(delineate_data_raw.temperature.data));
	xmiss = delineate_data_raw.transmissometer.data((delineate_data_raw.transmissometer.stationNo == stat) & ...
						        (delineate_data_raw.transmissometer.castNo == cast) & ...
						   	((delineate_data_raw.transmissometer.flag == 1) | (delineate_data_raw.transmissometer.flag == 2)) & ...
						        (delineate_data_raw.transmissometer.depth <= ppz) & ...
						   	~isnan(delineate_data_raw.transmissometer.data));
	po4 = delineate_data_raw.po4.data((int64(delineate_data_raw.po4.stationNo) == stat) & ...
					  (delineate_data_raw.po4.castNo == cast) & ...
					  ((delineate_data_raw.po4.flag == 1) | (delineate_data_raw.po4.flag == 2)) & ...
					  (delineate_data_raw.po4.depth <= ppz) & ...
					  ~isnan(delineate_data_raw.po4.data));
	no3 = delineate_data_raw.no3.data((int64(delineate_data_raw.no3.stationNo) == stat) & ...
					  (delineate_data_raw.no3.castNo == cast) & ...
					  ((delineate_data_raw.no3.flag == 1) | (delineate_data_raw.po4.flag == 2)) & ...
					  (delineate_data_raw.no3.depth <= ppz) & ...
					  ~isnan(delineate_data_raw.no3.data));
	pn = delineate_data_raw.pn.data((delineate_data_raw.pn.stationNo == stat) & ...
					  (delineate_data_raw.pn.depth <= ppz) & ...
					  ~isnan(delineate_data_raw.pn.data));

	% get mean of station ::
	temperature(istat) = mean(temp); 
	transmissometer(istat) = mean(xmiss); 
	phosphate(istat) = nanmean(po4); 
	nitrate(istat) = nanmean(no3); 
	particulate_n(istat) = mean(pn); 

end

%%  make data array ::
%   make initial ::
delineate_data.stationNo = gp15_stations.stationNo; 
delineate_data.latitude = gp15_stations.latitude; 
delineate_data.longitude = gp15_stations.longitude; 
delineate_data.temperature = temperature; 
delineate_data.transmissometer = transmissometer; 
delineate_data.phosphate = phosphate; 
delineate_data.nitrate = nitrate; 
delineate_data.particulate_n = particulate_n; 

%   find hplc indices and add in ::
[~, idxHplc, idxGp15] = intersect(gp15_hplc.station, gp15_stations.stationNo); 
delineate_data.micro_hplc = NaN(size(gp15_stations.stationNo)); 
delineate_data.nano_hplc = NaN(size(gp15_stations.stationNo)); 
delineate_data.pico_hplc = NaN(size(gp15_stations.stationNo)); 
delineate_data.micro_hplc(idxGp15) = gp15_hplc.micro;
delineate_data.nano_hplc(idxGp15) = gp15_hplc.nano;
delineate_data.pico_hplc(idxGp15) = gp15_hplc.pico;

%   make a table ::
delineate_data = struct2table(delineate_data); 

%%  plot
%   read in regions (from first looking at data) ::
regions = readtable([pro_output_basepath 'delineateRegions/regions/regions.xlsx']); 
k = (regions.first(2:end) + regions.last(1:end-1)) / 2; 

%   start plot ::
delineate_plot = figure; 
tiledlayout(3, 2, 'tileSpacing', 'compact'); 

%   temperature ::
nexttile(); 
scatter(gp15_stations.stationNo, temperature, 50, 'k', 'filled');
xline(k, 'linewidth', 1, 'linestyle', ':', 'color', 'black');
ylabel('Temp [$^\circ$C]', 'interpreter', 'latex');
title('(A) Temperature', 'interpreter', 'latex')
set(gca, 'xDir', 'reverse', 'box', 'on', 'tickLabelInterpreter', 'latex');

%   nitrate ::
nexttile(); 
scatter(gp15_stations.stationNo, nitrate, 50, 'k', 'filled');
xline(k, 'linewidth', 1, 'linestyle', ':', 'color', 'black');
ylabel('NO$_3$ [$\mu$M]', 'interpreter', 'latex');
title('(B) Nitrate', 'interpreter', 'latex')
set(gca, 'xDir', 'reverse', 'box', 'on', 'tickLabelInterpreter', 'latex');

%   phosphate ::
nexttile(); 
scatter(gp15_stations.stationNo, phosphate, 50, 'k', 'filled'); 
xline(k, 'linewidth', 1, 'linestyle', ':', 'color', 'black');
ylabel('PO$_4$ [$\mu$M]', 'interpreter', 'latex');
title('(C) Phosphate', 'interpreter', 'latex')
set(gca, 'xDir', 'reverse', 'box', 'on', 'tickLabelInterpreter', 'latex');

%   particulate nitrogen ::
nexttile(); 
scatter(gp15_stations.stationNo, particulate_n, 50, 'k', 'filled'); 
xline(k, 'linewidth', 1, 'linestyle', ':', 'color', 'black');
ylabel('P [$\mu$M]', 'interpreter', 'latex');
title('(D) Particulate N', 'interpreter', 'latex')
set(gca, 'xDir', 'reverse', 'box', 'on', 'tickLabelInterpreter', 'latex');

%   transmissometer ::
nexttile(); 
scatter(gp15_stations.stationNo, transmissometer, 50, 'k', 'filled');
xline(k, 'linewidth', 1, 'linestyle', ':', 'color', 'black');
ylabel({'Transmissometer', '[volts flipped]'}, 'interpreter', 'latex');
xlabel('Stations', 'interpreter', 'latex'); 
title('(E) Bulk Particles (Transm. Voltage)', 'interpreter', 'latex')
set(gca, 'xDir', 'reverse', 'yDir', 'reverse', 'box', 'on', 'tickLabelInterpreter', 'latex');

%   pigment ::
nexttile(); 
bar(gp15_hplc.station, 1E2 .* [gp15_hplc.micro'; gp15_hplc.nano'; gp15_hplc.pico']', 'stacked');
xline(k, 'linewidth', 1, 'linestyle', ':', 'color', 'black');
ylabel({'Pigment', 'Composition [\%]'}, 'interpreter', 'latex'); 
xlabel('Stations', 'interpreter', 'latex'); 
title('(F) Pigment Composition', 'interpreter', 'latex'); 
legend('Micro', 'Nano', 'Pico', 'location', 'eastOutside', 'interpreter', 'latex');
set(gca, 'xLim', [0, 40], 'xDir', 'reverse', 'yLim', [0, 1E2], 'box', 'on', 'tickLabelInterpreter', 'latex');

%   save ::
exportgraphics(delineate_plot, [plot_output_basepath 'delineateRegions/regionalMetrics.png'], 'resolution', 300);

%   save data ::
%%% .mat ::
save([pro_output_basepath 'delineateRegions/regionalData/delineate_data.mat'], 'delineate_data');
save([pro_output_basepath 'delineateRegions/regionalData/regions.mat'], 'regions');

%%% .xlsx ::
writetable(delineate_data, [pro_output_basepath 'delineateRegions/regionalData/delineate_data.xlsx']);

%%  end subroutine
