%%  load data
%   gp15 observations ::
load([pro_output_basepath 'doQc/gp15/gp15_obs.mat'], 'gp15_obs');

%   bottles data ::
load([pro_output_basepath 'readData/gp15/gp15_bottles.mat'], 'gp15_bottles');

%   ctd data ::
load([pro_output_basepath 'doQc/gp15/gp15_ctd.mat'], 'gp15_ctd');

%   ppz data ::
load([pro_output_basepath 'doQc/ppz/gp15_ppz.mat'], 'gp15_ppz');

%   station data:
load([pro_output_basepath 'doQc/stations/gp15_stations.mat'], 'gp15_stations', 'NUMSTAT');

%   mld data ::
load([sim_output_basepath 'calcMld/gp15_mld.mat'], 'gp15_mld');

%%  end subroutine
