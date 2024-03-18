%%  load data
%   inputs ::
load([pro_output_basepath 'collateData/gp15_inputs.mat'], 'gp15_inputs');

%   station data:
load([pro_output_basepath 'doQc/stations/gp15_stations.mat'], 'gp15_stations', 'NUMSTAT');

%   mld data ::
load([sim_output_basepath 'calcMld/gp15_mld.mat'], 'gp15_mld');

%   npp ::
load([pro_output_basepath 'doQc/gp15/gp15_npp.mat'], 'gp15_npp');

%   regions ::
load([pro_output_basepath 'delineateRegions/regionalData/regions.mat'], 'regions');
regionNames = {'North Pacific High Productivity Zone', 'North Pacific Gyre', 'Equatorial Pacific', 'South Pacific Gyre'}; 
numRegions = length(regionNames); 

%   ratio ::
load([sim_output_basepath 'regressRegionalRatio/ratioPocLsf.mat'], 'ratioPocLsf'); 
load([sim_output_basepath 'regressRegionalRatio/ratioPnLsf.mat'], 'ratioPnLsf'); 
load([sim_output_basepath 'regressRegionalRatio/ratioPocSsf.mat'], 'ratioPocSsf'); 
load([sim_output_basepath 'regressRegionalRatio/ratioPnSsf.mat'], 'ratioPnSsf'); 

%%  end subroutine
