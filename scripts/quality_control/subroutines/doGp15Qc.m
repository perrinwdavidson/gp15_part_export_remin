%%  load data
load([pro_output_basepath 'readData/gp15/gp15_obs.mat'], 'gp15_obs');
load([pro_output_basepath 'readData/gp15/gp15_part.mat'], 'gp15_part');
load([pro_output_basepath 'readData/gp15/gp15_ctd.mat'], 'gp15_ctd');
% load([pro_output_basepath 'readData/gp15/gp15_npp.mat'], 'gp15_npp');
load([pro_output_basepath 'readData/gp15/gp15_particles.mat'], 'gp15_particles');

%%  remove nans
gp15_obs = gp15_obs(~isnan(gp15_obs.th234), :);

%%  sort
gp15_obs = sortrows(gp15_obs, 'stationNo');

%%  only integer stations
gp15_obs = gp15_obs(mod(gp15_obs.stationNo, 1) == 0, :);

%%  make potential density anomaly
gp15_ctd.ctd_potentialDensity = gsw_sigma0(gp15_ctd.ctd_salinity, gp15_ctd.ctd_temperature);

%%  only keep ppz stations
idx_keep = ismember(gp15_obs.stationNo, gp15_stations.stationNo);
gp15_obsNew = gp15_obs(idx_keep, :); 

%%  generate station times 
statTimeCount;

%%  quality control the particulate data
idx_particles_stats = ismember(gp15_particles.station_num, gp15_stations.stationNo);
gp15_particles = gp15_particles(idx_particles_stats, :); 
gp15_particles.POC_LPT_uM(gp15_particles.POC_LPT_uM_flag ~= 1) = NaN;
gp15_particles.PN_LPT_uM(gp15_particles.PN_LPT_uM_flag ~= 1) = NaN;
gp15_particles.POC_SPT_uM(gp15_particles.POC_SPT_uM_flag ~= 1) = NaN;
gp15_particles.PN_SPT_uM(gp15_particles.PN_SPT_uM_flag ~= 1) = NaN;
gp15_particles = gp15_particles(:, [2 3 7:10 15:18]);

%%  average duplicate data
aveDuplicate; 

%%  quality control particulate data
qcGp15Poc;  
qcGp15Pn;

%%  save data
%   excel ::
writetable(gp15_obs, [pro_output_basepath 'doQc/gp15/gp15_obs.xlsx'], 'writeMode', 'overwritesheet');
writetable(gp15_obsNoQc, [pro_output_basepath 'doQc/gp15/gp15_obsNoQc.xlsx'], 'writeMode', 'overwritesheet');
writetable(gp15_part, [pro_output_basepath 'doQc/gp15/gp15_part.xlsx'], 'writeMode', 'overwritesheet');
writetable(gp15_ctd, [pro_output_basepath 'doQc/gp15/gp15_ctd.xlsx'], 'writeMode', 'overwritesheet');
writetable(gp15_particles, [pro_output_basepath 'doQc/gp15/gp15_particles.xlsx'], 'writeMode', 'overwritesheet');
writetable(gp15_stations, [pro_output_basepath 'doQc/stations/gp15_stations.xlsx'], 'writeMode', 'overwritesheet');

%   mat ::
save([pro_output_basepath 'doQc/gp15/gp15_obs.mat'], 'gp15_obs');
save([pro_output_basepath 'doQc/gp15/gp15_obsNoQc.mat'], 'gp15_obsNoQc');
save([pro_output_basepath 'doQc/gp15/gp15_part.mat'], 'gp15_part');
save([pro_output_basepath 'doQc/gp15/gp15_ctd.mat'], 'gp15_ctd');
save([pro_output_basepath 'doQc/gp15/gp15_particles.mat'], 'gp15_particles');
save([pro_output_basepath 'doQc/stations/gp15_stations.mat'], 'gp15_stations', 'NUMSTAT');

%% end subroutine
