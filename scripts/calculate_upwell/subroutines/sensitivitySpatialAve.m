%%  sensitivity: spatial averaging window for upwelling velocity
%   tests four spatial averaging levels — no averaging (window = 1) and
%   three progressively wider windows including the nominal value — with
%   temporal averaging fixed at deltaDay = 35 d. runs inside the calcUpwell
%   product loop (and the ECCO block) after spatialTimeAverage. ::

%%  apply temporal averaging only (holds temporal constant while varying spatial)
wTimeOnly = movmean(w, timeAve, 4, 'omitnan');

%%  define spatial averaging test levels
%   [1 = none, ~half nominal, nominal, ~double nominal]; clamped to int >= 1.
%   unique() removes duplicates that arise for small spatAve (e.g. ECCO). ::
spatAveTest = unique(max(1, round([1, spatAve / 2, spatAve, spatAve * 2])));
nTest       = length(spatAveTest);
lineStyles  = {'-', '--', ':', '-.'};

%%  apply each spatial window to the temporally-averaged field
%   window = 1 is a true identity (movmean with window 1 is a no-op), but
%   bypassing movmean entirely makes the no-averaging intent unambiguous. ::
wSens = cell(nTest, 1);
for iSp = 1 : 1 : nTest
    if spatAveTest(iSp) == 1
        wSens{iSp} = wTimeOnly;
    else
        wTmp       = movmean(wTimeOnly, spatAveTest(iSp), 1, 'omitnan');
        wSens{iSp} = movmean(wTmp,     spatAveTest(iSp), 2, 'omitnan');
    end
end

%%  extract comparison slice
%   100 m depth, time nearest mean cruise date, mean GP15 station longitude ::
idxDepth100   = find(abs(u_mercator.depth - 100) == min(abs(u_mercator.depth - 100)), 1);
idxTimeCruise = find(abs(u_mercator.time - mean(gp15_stations.date)) == min(abs(u_mercator.time - mean(gp15_stations.date))), 1);
idxLonMean    = find(abs(u_mercator.longitude - mean(mod(gp15_stations.longitude + 360, 360))) == min(abs(u_mercator.longitude - mean(mod(gp15_stations.longitude + 360, 360)))), 1);

%%  plot
figure;
hold('on');
for iSp = 1 : 1 : nTest
    wSlice = squeeze(wSens{iSp}(idxLonMean, :, idxDepth100, idxTimeCruise)) * SEC2DAY;
    if spatAveTest(iSp) == 1
        lbl = 'No spatial avg';
    else
        lbl = ['$' num2str(spatAveTest(iSp) * spaceResolution) '^\circ$ spatial avg'];
    end
    plot(u_mercator.latitude, wSlice, lineStyles{iSp}, 'lineWidth', 2, 'displayName', lbl);
end
hold('off');
xline(0, '-k', 'lineWidth', 1, 'handleVisibility', 'off');
xlabel('\textbf{Latitude [deg N.]}', 'interpreter', 'latex');
ylabel('$w$ \textbf{[m d$^{-1}$]}', 'interpreter', 'latex');
legend('location', 'best', 'interpreter', 'latex', 'fontSize', 14);
set(gca, 'tickLabelInterpreter', 'latex', 'fontSize', 16, 'fontWeight', 'bold', 'lineWidth', 1, 'box', 'on');
title(['\textbf{Sensitivity: Spatial Averaging Window (' upper(dataProduct) ', 100 m, 35 d)}'], ...
      'interpreter', 'latex', 'fontSize', 20);
set(gcf, 'position', [0, 0, 1200, 500]);
exportgraphics(gcf, [plot_output_basepath 'calcUpwell/sensitivity_spatAve_' dataProduct '.pdf'], ...
               'ContentType', 'vector');

%%  end subroutine
