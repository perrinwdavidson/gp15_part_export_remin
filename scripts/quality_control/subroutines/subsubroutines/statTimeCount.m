%% sample time and count for stations
for iTime = 1 : 1 : NUMSTAT
    
    sn = gp15_stations.stationNo(iTime);
    dat = gp15_obs.date(gp15_obs.stationNo == sn);
    gp15_stations.date(iTime) = dateshift(mean(dat), 'start', 'day', 'nearest');     

end

%% end subsubroutine
