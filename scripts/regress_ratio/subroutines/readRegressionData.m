%% read all data
%  observations ::
load([pro_output_basepath 'interpData/observations/gp15_obs.mat'], 'gp15_obs');
gp15_obsFull = gp15_obs; 
clear('gp15_obs'); 
load([pro_output_basepath 'doQc/gp15/gp15_obs.mat'], 'gp15_obs');

%   ppz data ::
load([pro_output_basepath 'doQc/ppz/gp15_ppz.mat'], 'gp15_ppz');

%   mld data ::
load([sim_output_basepath 'calcMld/gp15_mld.mat'], 'gp15_mld');

%   station data ::
load([pro_output_basepath 'doQc/stations/gp15_stations.mat'], 'gp15_stations', 'NUMSTAT');

%  regions ::
load([pro_output_basepath 'delineateRegions/regionalData/regions.mat'], 'regions');

%  set region names ::
region_names = {'North Pacific High Productivity Zone', 'North Pacific Gyre', 'Equatorial Pacific', 'South Pacific Gyre'}; 

%% end subroutine
