%%  load data
%   gp15 observations ::
load([pro_output_basepath 'doQc/gp15/gp15_obs.mat'], 'gp15_obs');

%   ppz data ::
load([pro_output_basepath 'doQc/ppz/gp15_ppz.mat'], 'gp15_ppz');

%   mld data ::
load([sim_output_basepath 'calcMld/gp15_mld.mat'], 'gp15_mld');

%   eqp data ::
load([sim_output_basepath 'calcEqp/gp15_eqp.mat'], 'gp15_eqp');

%   station data:
load([pro_output_basepath 'doQc/stations/gp15_stations.mat'], 'gp15_stations', 'NUMSTAT');

%%  end subroutine
