%%  load data 
%   load gp15 data ::
load([pro_output_basepath 'doQc/stations/gp15_stations.mat'], 'gp15_stations', 'NUMSTAT');

%   load npp data ::
load([pro_output_basepath 'readData/cbpm/npp.mat'], 'NPP');

%   load npp data ::
X = ncread([sim_output_basepath 'calcNpp/decorrelation/npp_matrix.nc'], 'longitude');
Y = ncread([sim_output_basepath 'calcNpp/decorrelation/npp_matrix.nc'], 'latitude');
T = ncread([sim_output_basepath 'calcNpp/decorrelation/npp_matrix.nc'], 'time');
V = ncread([sim_output_basepath 'calcNpp/decorrelation/npp_matrix.nc'], 'npp');

%   load decorrelation scales ::
lxRaw = readmatrix([sim_output_basepath 'calcNpp/decorrelation/spatial_decorrelation.csv']);
ltRaw = readmatrix([sim_output_basepath 'calcNpp/decorrelation/temporal_decorrelation.csv']);

%   appropriate round ::
lx = km2deg(round(lxRaw(2, 1),  1 + int64(floor(log10(lxRaw(2, 2)))), 'significant'));  % [deg]
lt = round(ltRaw(2, 1),  1 + int64(floor(log10(ltRaw(2, 2)))), 'significant');  % [d]

%%  end subroutine
