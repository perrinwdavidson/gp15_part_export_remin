%%  extract w slice from spatially and temporally averaged field
%   indices: 100 m depth, mean GP15 cruise date, mean GP15 transect longitude.
%   uses wSpatAve (35-d temporal mean applied in spatialTimeAverage.m) —
%   the field that enters the flux model — at a single depth level.
%   old version used raw w + depth-mean + one movmean pass, which produced
%   a noisy result; those three lines are commented below for reference.
% idxDepth = find(min(abs(u_mercator.depth - 100)) == abs(u_mercator.depth - 100));
% wAveDay  = squeeze(mean(w(idxLon, :, 1:idxDepth, idxTime), 3))' * SEC2DAY;
% wAveDay  = movmean(wAveDay, degreeAve / spaceResolution, 1);
%--------------------------------------------------------------------------
idxLon   = find(abs(u_mercator.longitude - mean(mod(gp15_stations.longitude + 360, 360))) == min(abs(u_mercator.longitude - mean(mod(gp15_stations.longitude + 360, 360)))), 1);
idxDepth = find(abs(u_mercator.depth - 100) == min(abs(u_mercator.depth - 100)), 1);
idxTime  = find(abs(u_mercator.time - mean(gp15_stations.date)) == min(abs(u_mercator.time - mean(gp15_stations.date))), 1);

wAveDay = squeeze(wSpatAve(idxLon, :, idxDepth, idxTime)) * SEC2DAY;

%%  figure 1: meridional slice at mean GP15 longitude and mean cruise date
figure;
plot(u_mercator.latitude, wAveDay, '-k', 'lineWidth', 1.5);
xline(0, '-k', 'lineWidth', 0.5, 'handleVisibility', 'off');
xlabel('\textbf{Latitude ($^{\circ}$N)}', 'interpreter', 'latex', 'fontSize', 18);
ylabel('\textbf{Vertical Velocity (m day$^{-1}$)}', 'interpreter', 'latex', 'fontSize', 18);
title(['\textbf{Spatially and Temporally Averaged Upwelling at 100 m (' upper(dataProduct) ', $' num2str(spatAve * spaceResolution) '^\circ$ spatial avg, 35 d)}'], ...
      'interpreter', 'latex', 'fontSize', 18);
set(gca, 'tickLabelInterpreter', 'latex', 'fontSize', 16, 'fontWeight', 'bold', 'lineWidth', 1, 'box', 'on');
set(gcf, 'units', 'inches', 'position', [0, 0, 20, 8], 'paperUnits', 'inches', 'paperSize', [20, 8]);
exportgraphics(gcf, [plot_output_basepath 'calcUpwell/w_' dataProduct '.pdf'], 'ContentType', 'vector');

%%  figure 2: wSpatAve sampled at closest grid point to each GP15 station
%   for each station, finds the nearest (lon, lat, time) index in wSpatAve
%   and extracts w at 100 m. this is the model value at the actual station
%   position and sampling date — directly comparable to the kriging output
%   from interpUpwell.m, which queries the same coordinates. ::
wStation = NaN(NUMSTAT, 1);
for iStat = 1 : 1 : NUMSTAT
    idxLon_i  = find(abs(u_mercator.longitude - mod(gp15_stations.longitude(iStat) + 360, 360)) == ...
                     min(abs(u_mercator.longitude - mod(gp15_stations.longitude(iStat) + 360, 360))), 1);
    idxLat_i  = find(abs(u_mercator.latitude  - gp15_stations.latitude(iStat)) == ...
                     min(abs(u_mercator.latitude  - gp15_stations.latitude(iStat))), 1);
    idxTime_i = find(abs(u_mercator.time - gp15_stations.date(iStat)) == ...
                     min(abs(u_mercator.time - gp15_stations.date(iStat))), 1);
    wStation(iStat) = wSpatAve(idxLon_i, idxLat_i, idxDepth, idxTime_i) * SEC2DAY;
end

figure;
plot(gp15_stations.latitude, wStation, '-ok', 'lineWidth', 1.5, ...
     'markerFaceColor', 'k', 'markerSize', 10);
xline(0, '-k', 'lineWidth', 0.5, 'handleVisibility', 'off');
xlabel('\textbf{Latitude ($^{\circ}$N)}', 'interpreter', 'latex', 'fontSize', 18);
ylabel('\textbf{Vertical Velocity (m day$^{-1}$)}', 'interpreter', 'latex', 'fontSize', 18);
title(['\textbf{$w$ at GP15 Station Positions (' upper(dataProduct) ', 100 m, $' num2str(spatAve * spaceResolution) '^\circ$ spatial avg, 35 d)}'], ...
      'interpreter', 'latex', 'fontSize', 18);
set(gca, 'tickLabelInterpreter', 'latex', 'fontSize', 16, 'fontWeight', 'bold', 'lineWidth', 1, 'box', 'on');
set(gcf, 'units', 'inches', 'position', [0, 0, 20, 8], 'paperUnits', 'inches', 'paperSize', [20, 8]);
exportgraphics(gcf, [plot_output_basepath 'calcUpwell/w_' dataProduct '_stations.pdf'], ...
               'ContentType', 'vector');

%%  end subroutine
