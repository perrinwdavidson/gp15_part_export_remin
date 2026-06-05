%%  read data
%   ctd data ::
load([pro_output_basepath 'doQc/gp15/gp15_ctd.mat'], 'gp15_ctd');

%   bottle data ::
load([pro_output_basepath 'readData/gp15/gp15_bottles.mat'], 'gp15_bottles');

%   observations ::
load([pro_output_basepath 'interpData/observations/gp15_obs.mat'], 'gp15_obs');

%   hplc ::
load([pro_output_basepath 'readData/gp15/gp15_hplc.mat'], 'gp15_hplc');

%   ppz data ::
load([pro_output_basepath 'doQc/ppz/gp15_ppz.mat'], 'gp15_ppz');

%   mld data ::
load([sim_output_basepath 'calcMld/gp15_mld.mat'], 'gp15_mld');

%   station data ::
load([pro_output_basepath 'doQc/stations/gp15_stations.mat'], 'gp15_stations', 'NUMSTAT');

%%  end subroutine
