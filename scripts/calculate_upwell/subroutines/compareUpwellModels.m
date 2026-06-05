%%  compareUpwellModels - compare CMEMS ensemble members, their mean, and GREP
%   runs inside calcUpwell.m after the dataProducts loop and ECCO block.
%   loads the saved w netCDFs for cglo, foam, glor, oras, and grep.
%   the saved w is already wSpatAve (spatially + 35-d temporally averaged),
%   so no further smoothing is applied after extraction.
%   produces two figures:
%     figure 1: meridional slice at mean GP15 longitude and mean cruise date
%     figure 2: wSpatAve sampled at closest grid point to each GP15 station
%               (lon_i, lat_i, 100 m, date_i) — directly comparable to
%               interpUpwell.m output.
%--------------------------------------------------------------------------

%%  define members
members      = {'cglo', 'foam', 'glor', 'oras'};
memberLabels = {'CGLO', 'FOAM', 'GLOR', 'ORAS'};
nMembers     = length(members);

%%  CMEMS spatial averaging window (spaceResolution was overwritten to 1.0
%   by the ECCO block; recompute from the known 0.25-deg CMEMS grid) ::
cmems_spaceResolution = 0.25;
cmems_spatAve = floor(degreeAve / cmems_spaceResolution);

%%  load CMEMS grid coordinates from saved netCDF
%   u_mercator is overwritten by the ECCO block before compareUpwellModels
%   runs, so its longitude/latitude/time now reflect the ECCO 1-deg global
%   grid — not the CMEMS 0.25-deg regional grid. load coordinates directly
%   from w_grep.nc (all CMEMS products share the same grid). time is stored
%   as datenum in the netCDF; convert back to datetime for comparison. ::
refFile    = [sim_output_basepath 'calcUpwell/w_grep.nc'];
cmems_lon  = ncread(refFile, 'longitude');
cmems_lat  = ncread(refFile, 'latitude');
cmems_dep  = ncread(refFile, 'depth');
cmems_time = datetime(ncread(refFile, 'time'), 'convertFrom', 'datenum');

%%  shared grid indices (all CMEMS products share the same 0.25 deg grid)

%   mean-slice indices ::
idxLon   = find(abs(cmems_lon  - mean(mod(gp15_stations.longitude + 360, 360))) == min(abs(cmems_lon  - mean(mod(gp15_stations.longitude + 360, 360)))), 1);
idxDepth = find(abs(cmems_dep  - 100)                                           == min(abs(cmems_dep  - 100)),                                           1);
idxTime  = find(abs(cmems_time - mean(gp15_stations.date))                      == min(abs(cmems_time - mean(gp15_stations.date))),                      1);
wLat     = cmems_lat;

%   per-station indices (one set of lon/lat/time indices per GP15 station) ::
idxLon_s  = NaN(NUMSTAT, 1);
idxLat_s  = NaN(NUMSTAT, 1);
idxTime_s = NaN(NUMSTAT, 1);
for iStat = 1 : 1 : NUMSTAT
    idxLon_s(iStat)  = find(abs(cmems_lon  - mod(gp15_stations.longitude(iStat) + 360, 360)) == min(abs(cmems_lon  - mod(gp15_stations.longitude(iStat) + 360, 360))), 1);
    idxLat_s(iStat)  = find(abs(cmems_lat  - gp15_stations.latitude(iStat))                  == min(abs(cmems_lat  - gp15_stations.latitude(iStat))),                  1);
    idxTime_s(iStat) = find(abs(cmems_time - gp15_stations.date(iStat))                       == min(abs(cmems_time - gp15_stations.date(iStat))),                       1);
end

%%  load each member once; extract both mean-slice and per-station values
wMembers        = NaN(length(wLat), nMembers);
wMembersStation = NaN(NUMSTAT,      nMembers);
for iMember = 1 : 1 : nMembers
    wRaw = ncread([sim_output_basepath 'calcUpwell/w_' members{iMember} '.nc'], 'w');
    wMembers(:, iMember) = squeeze(wRaw(idxLon, :, idxDepth, idxTime)) * SEC2DAY;
    for iStat = 1 : 1 : NUMSTAT
        wMembersStation(iStat, iMember) = wRaw(idxLon_s(iStat), idxLat_s(iStat), idxDepth, idxTime_s(iStat)) * SEC2DAY;
    end
end

%%  ensemble means
wEnsembleMean        = mean(wMembers,        2, 'omitnan');
wEnsembleMeanStation = mean(wMembersStation, 2, 'omitnan');

%%  load GREP; extract both slices
wRawGrep     = ncread([sim_output_basepath 'calcUpwell/w_grep.nc'], 'w');
wGrep        = squeeze(wRawGrep(idxLon, :, idxDepth, idxTime)) * SEC2DAY;
wGrepStation = NaN(NUMSTAT, 1);
for iStat = 1 : 1 : NUMSTAT
    wGrepStation(iStat) = wRawGrep(idxLon_s(iStat), idxLat_s(iStat), idxDepth, idxTime_s(iStat)) * SEC2DAY;
end

%%  shared styling
memberLineStyles = {'-', '--', '-.', ':'};
memberColor      = [0.65, 0.65, 0.65];

%%  figure 1: meridional slice at mean GP15 longitude and mean cruise date
figure;
hold('on');
for iMember = 1 : 1 : nMembers
    plot(wLat, wMembers(:, iMember), memberLineStyles{iMember}, ...
         'color', memberColor, 'lineWidth', 1, 'displayName', memberLabels{iMember});
end
plot(wLat, wEnsembleMean, '-k',  'lineWidth', 2, 'displayName', 'CMEMS Mean');
plot(wLat, wGrep,         '--r', 'lineWidth', 2, 'displayName', 'GREP');
hold('off');
box('on');
xlabel('\textbf{Latitude ($^{\circ}$N)}', 'interpreter', 'latex', 'fontSize', 20);
ylabel('\textbf{Vertical Velocity (m day$^{-1}$)}', 'interpreter', 'latex', 'fontSize', 20);
legend('location', 'best', 'interpreter', 'latex', 'fontSize', 14);
title(['\textbf{CMEMS Ensemble Members, Mean, and GREP (Meridional Slice, 100 m, $' num2str(cmems_spatAve * cmems_spaceResolution) '^\circ$ spatial avg, 35 d)}'], ...
      'interpreter', 'latex', 'fontSize', 18);
set(gca, 'tickLabelInterpreter', 'latex', 'fontSize', 16, 'fontWeight', 'bold', 'lineWidth', 1, 'box', 'on');
set(gcf, 'units', 'inches', 'position', [0, 0, 20, 8], 'paperUnits', 'inches', 'paperSize', [20, 8]);
exportgraphics(gcf, [plot_output_basepath 'calcUpwell/w_comparison.pdf'], 'ContentType', 'vector');

%%  figure 2: wSpatAve at closest grid point to each GP15 station
figure;
hold('on');
for iMember = 1 : 1 : nMembers
    plot(gp15_stations.latitude, wMembersStation(:, iMember), memberLineStyles{iMember}, ...
         'color', memberColor, 'lineWidth', 1, 'displayName', memberLabels{iMember});
end
plot(gp15_stations.latitude, wEnsembleMeanStation, '-k',  'lineWidth', 2, 'displayName', 'CMEMS Mean');
plot(gp15_stations.latitude, wGrepStation,         '--r', 'lineWidth', 2, 'displayName', 'GREP');
hold('off');
box('on');
xlabel('\textbf{Latitude ($^{\circ}$N)}', 'interpreter', 'latex', 'fontSize', 20);
ylabel('\textbf{Vertical Velocity (m day$^{-1}$)}', 'interpreter', 'latex', 'fontSize', 20);
legend('location', 'best', 'interpreter', 'latex', 'fontSize', 14);
title(['\textbf{CMEMS Ensemble Members, Mean, and GREP (GP15 Stations, 100 m, $' num2str(cmems_spatAve * cmems_spaceResolution) '^\circ$ spatial avg, 35 d)}'], ...
      'interpreter', 'latex', 'fontSize', 18);
set(gca, 'tickLabelInterpreter', 'latex', 'fontSize', 16, 'fontWeight', 'bold', 'lineWidth', 1, 'box', 'on');
set(gcf, 'units', 'inches', 'position', [0, 0, 20, 8], 'paperUnits', 'inches', 'paperSize', [20, 8]);
exportgraphics(gcf, [plot_output_basepath 'calcUpwell/w_comparison_stations.pdf'], 'ContentType', 'vector');

%%  end subroutine
