%% plotQcPn - diagnostic 8-panel PN QC plots per station
%  standalone; run after doQc.m has saved gp15_obsNoQc.mat.
%  uses pre-flier-removal data (gp15_obsNoQc) to reproduce the original
%  QC figures, showing interpolated profiles and raw ratios before flier
%  correction.
%--------------------------------------------------------------------------
setupModel;

%%  load saved outputs
load([pro_output_basepath 'doQc/gp15/gp15_obs.mat'], 'gp15_obs');
load([pro_output_basepath 'doQc/stations/gp15_stations.mat'], 'gp15_stations', 'NUMSTAT');

%%  loop over stations
for iStat = 1 : 1 : NUMSTAT

    statNo   = gp15_stations.stationNo(iStat);
    statData = gp15_obs(gp15_obs.stationNo == statNo, :);

    figure;
    tl = tiledlayout(2, 4, 'tileSpacing', 'compact');

    % u238 ::
    nexttile();
    errorbar(statData.u238, statData.depth, [], [], statData.uncertU238, statData.uncertU238, '-ok', ...
             'markerEdgeColor', 'k', 'markerFaceColor', 'k', 'lineWidth', 1, 'markerSize', 7.5, ...
             'capSize', 0);
    box('on');
    xlabel('\textbf{$^{238}$U (dpm L$^{-1}$)}', 'interpreter', 'latex', 'fontSize', 20);
    ylabel('\textbf{Depth (m)}', 'interpreter', 'latex', 'fontSize', 20);
    ylim([0 1000]);
    set(gca, 'yDir', 'reverse', 'tickLabelInterpreter', 'latex', 'fontSize', 16, 'fontWeight', 'bold', 'lineWidth', 1);

    % th234 ::
    nexttile();
    errorbar(statData.th234, statData.depth, [], [], statData.uncertTh234, statData.uncertTh234, '-ok', ...
             'markerEdgeColor', 'k', 'markerFaceColor', 'k', 'lineWidth', 1, 'markerSize', 7.5, ...
             'capSize', 0);
    box('on');
    xlabel('\textbf{Total $^{234}$Th (dpm L$^{-1}$)}', 'interpreter', 'latex', 'fontSize', 20);
    ylim([0 1000]);
    set(gca, 'yDir', 'reverse', 'tickLabelInterpreter', 'latex', 'fontSize', 16, 'fontWeight', 'bold', 'lineWidth', 1);

    % th234 ssf ::
    nexttile();
    errorbar(statData.th234PocSmall, statData.depth, [], [], statData.uncertTh234PocSmall, statData.uncertTh234PocSmall, '-ok', ...
             'markerEdgeColor', 'k', 'markerFaceColor', 'k', 'lineWidth', 1, 'markerSize', 7.5, ...
             'capSize', 0);
    box('on');
    xlabel('\textbf{SSF Particulate $^{234}$Th (dpm L$^{-1}$)}', 'interpreter', 'latex', 'fontSize', 20);
    ylim([0 1000]);
    set(gca, 'yDir', 'reverse', 'tickLabelInterpreter', 'latex', 'fontSize', 16, 'fontWeight', 'bold', 'lineWidth', 1);

    % th234 lsf ::
    nexttile();
    errorbar(statData.th234PocLarge, statData.depth, [], [], statData.uncertTh234PocLarge, statData.uncertTh234PocLarge, '-ok', ...
             'markerEdgeColor', 'k', 'markerFaceColor', 'k', 'lineWidth', 1, 'markerSize', 7.5, ...
             'capSize', 0);
    box('on');
    xlabel('\textbf{LSF Particulate $^{234}$Th (dpm L$^{-1}$)}', 'interpreter', 'latex', 'fontSize', 20);
    ylim([0 1000]);
    set(gca, 'yDir', 'reverse', 'tickLabelInterpreter', 'latex', 'fontSize', 16, 'fontWeight', 'bold', 'lineWidth', 1);

    % ssf pn ::
    nexttile();
    errorbar(statData.pnSmall, statData.depth, [], [], statData.uncertPnSmall, statData.uncertPnSmall, '-ok', ...
             'markerEdgeColor', 'k', 'markerFaceColor', 'k', 'lineWidth', 1, 'markerSize', 7.5, ...
             'capSize', 0);
    box('on');
    xlabel('\textbf{SSF PN ($\mu$M)}', 'interpreter', 'latex', 'fontSize', 20);
    ylabel('\textbf{Depth (m)}', 'interpreter', 'latex', 'fontSize', 20);
    ylim([0 1000]);
    set(gca, 'yDir', 'reverse', 'tickLabelInterpreter', 'latex', 'fontSize', 16, 'fontWeight', 'bold', 'lineWidth', 1);

    % lsf pn ::
    nexttile();
    errorbar(statData.pnLarge, statData.depth, [], [], statData.uncertPnLarge, statData.uncertPnLarge, '-ok', ...
             'markerEdgeColor', 'k', 'markerFaceColor', 'k', 'lineWidth', 1, 'markerSize', 7.5, ...
             'capSize', 0);
    box('on');
    xlabel('\textbf{LSF PN ($\mu$M)}', 'interpreter', 'latex', 'fontSize', 20);
    ylim([0 1000]);
    set(gca, 'yDir', 'reverse', 'tickLabelInterpreter', 'latex', 'fontSize', 16, 'fontWeight', 'bold', 'lineWidth', 1);

    % ssf pn:234th ::
    nexttile();
    errorbar(statData.pnTh234RatioSmall, statData.depth, [], [], statData.uncertPnTh234RatioSmall, statData.uncertPnTh234RatioSmall, '-ok', ...
             'markerEdgeColor', 'k', 'markerFaceColor', 'k', 'lineWidth', 1, 'markerSize', 7.5, ...
             'capSize', 0);
    box('on');
    xlabel('\textbf{SSF PN:$^{234}$Th ($\mu$mol dpm$^{-1}$)}', 'interpreter', 'latex', 'fontSize', 20);
    ylim([0 1000]);
    set(gca, 'yDir', 'reverse', 'tickLabelInterpreter', 'latex', 'fontSize', 16, 'fontWeight', 'bold', 'lineWidth', 1);

    % lsf pn:234th ::
    nexttile();
    errorbar(statData.pnTh234RatioLarge, statData.depth, [], [], statData.uncertPnTh234RatioLarge, statData.uncertPnTh234RatioLarge, '-ok', ...
             'markerEdgeColor', 'k', 'markerFaceColor', 'k', 'lineWidth', 1, 'markerSize', 7.5, ...
             'capSize', 0);
    box('on');
    xlabel('\textbf{LSF PN:$^{234}$Th ($\mu$mol dpm$^{-1}$)}', 'interpreter', 'latex', 'fontSize', 20);
    ylim([0 1000]);
    set(gca, 'yDir', 'reverse', 'tickLabelInterpreter', 'latex', 'fontSize', 16, 'fontWeight', 'bold', 'lineWidth', 1);

    title(tl, ['\textbf{Station ' num2str(statNo) '}'], 'interpreter', 'latex', 'fontSize', 24);
    set(gcf, 'units', 'inches', 'position', [0, 0, 20, 15], 'paperUnits', 'inches', 'paperSize', [20, 15]);
    exportgraphics(gcf, [plot_output_basepath, 'doQc/stationPlots/pn/station', num2str(statNo), '.pdf'], ...
                   'ContentType', 'vector');

end
close('all');

%%  end program
