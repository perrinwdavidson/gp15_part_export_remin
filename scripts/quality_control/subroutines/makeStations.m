%%  make array
stations = [gp15_ppz.stationNo, gp15_ppz.castNo, gp15_ppz.longitude, ...
            gp15_ppz.latitude, gp15_ppz.ppzDepth, gp15_ppz.normalizationTo100Meters, gp15_ppz.maxDepth400Meters];
stations = array2table(stations);
stations.Properties.VariableNames = {'stationNo', 'castNo', 'longitude', ...
                                     'latitude', 'depthPpz', 'depth100', 'depth400'};

%%  number of stations
NUMSTAT = size(stations, 1);

%%  save
%   give new variables ::
gp15_stations = stations;
clear('stations');

%   excel ::
writetable(gp15_stations, [pro_output_basepath 'doQc/stations/gp15_stations.xlsx']);

%   mat ::
save([pro_output_basepath 'doQc/stations/gp15_stations.mat'], 'gp15_stations', 'NUMSTAT');

%% end subroutine
