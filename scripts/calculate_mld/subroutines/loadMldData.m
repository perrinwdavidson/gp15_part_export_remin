%%  read data
%   ctd data ::
load([pro_output_basepath 'doQc/gp15/gp15_ctd.mat'], 'gp15_ctd');
load([pro_output_basepath 'readData/gp15/gp15_bottles.mat'], 'gp15_bottles');

%   ppz data ::
load([pro_output_basepath 'doQc/ppz/gp15_ppz.mat'], 'gp15_ppz');

%   station data ::
load([pro_output_basepath 'doQc/stations/gp15_stations.mat'], 'gp15_stations', 'NUMSTAT');

%   mld guesses ::
fileName = [input_basepath 'mld/mld_calcs_jak.xlsx'];
mldCalcs = readtable(fileName, 'variableNamingRule', 'preserve');

%%  end subroutine
