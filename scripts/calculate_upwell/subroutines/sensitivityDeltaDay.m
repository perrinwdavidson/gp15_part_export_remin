%%  sensitivity: temporal averaging window for upwelling velocity
%   computes wSpatAve for deltaDays = [10, 24, 35] d using the same spatial
%   averaging already applied in spatialTimeAverage.m.  runs inside the
%   calcUpwell product loop immediately after spatialTimeAverage. ::

%%  compute spatially-averaged field (no temporal averaging applied yet)
%   replicates the two spatial passes from spatialTimeAverage.m ::
wSpatOnly = movmean(w, spatAve, 1, 'omitnan');
wSpatOnly = movmean(wSpatOnly, spatAve, 2, 'omitnan');

%   temporal windows to test [d] ::
deltaDaysTest = [10, 24, 35];
nTest         = length(deltaDaysTest);

%   apply each temporal window ::
wSens = cell(nTest, 1);
for iDelta = 1 : 1 : nTest
    timeAveTest   = floor(deltaDaysTest(iDelta) / timeResolution);
    wSens{iDelta} = movmean(wSpatOnly, timeAveTest, 4, 'omitnan');
end

%%  extract comparison slice
%   100 m depth, time nearest mean cruise date, mean GP15 station longitude ::
idxDepth100   = find(abs(u_mercator.depth - 100)                                        == min(abs(u_mercator.depth - 100)),                                        1);
idxTimeCruise = find(abs(u_mercator.time  - mean(gp15_stations.date))                   == min(abs(u_mercator.time  - mean(gp15_stations.date))),                   1);
idxLonMean    = find(abs(u_mercator.longitude - mean(mod(gp15_stations.longitude + 360, 360))) == min(abs(u_mercator.longitude - mean(mod(gp15_stations.longitude + 360, 360)))), 1);

%%  plot
figure;
hold('on');
lineStyles = {'-', '--', ':'};
for iDelta = 1 : 1 : nTest
    wSlice = squeeze(wSens{iDelta}(idxLonMean, :, idxDepth100, idxTimeCruise)) * SEC2DAY;
    plot(u_mercator.latitude, wSlice, lineStyles{iDelta}, 'lineWidth', 2, ...
         'displayName', [num2str(deltaDaysTest(iDelta)) ' d']);
end
hold('off');
xline(0, '-k', 'lineWidth', 1, 'handleVisibility', 'off');
xlabel('\textbf{Latitude [deg N.]}', 'interpreter', 'latex');
ylabel('$w$ \textbf{[m d$^{-1}$]}', 'interpreter', 'latex');
legend('location', 'best', 'interpreter', 'latex', 'fontSize', 14);
set(gca, 'tickLabelInterpreter', 'latex', 'fontSize', 16, 'fontWeight', 'bold', 'lineWidth', 1, 'box', 'on');
title(['\textbf{Sensitivity: Temporal Averaging Window (' upper(dataProduct) ', 100 m, $' num2str(spatAve * spaceResolution) '^\circ$ spatial avg)}'], ...
      'interpreter', 'latex', 'fontSize', 20);
set(gcf, 'position', [0, 0, 1200, 500]);
exportgraphics(gcf, [plot_output_basepath 'calcUpwell/sensitivity_deltaDay_' dataProduct '.pdf'], 'ContentType', 'vector');

%%  end subroutine
