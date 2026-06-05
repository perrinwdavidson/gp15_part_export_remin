%% plotEqp - diagnostic plots of equilibrium point depth profiles
%  standalone; run after calcEqp.m has saved gp15_eqp.mat.
%  plots Th-234 and U-238 activity profiles with EQP, MLD, and PPZ markers
%  per station.
%--------------------------------------------------------------------------
setupModel;
setModelCoefficients;
configureCalcEqp;

%%  load saved outputs
load([sim_output_basepath 'calcEqp/gp15_eqp.mat'],          'gp15_eqp');
load([pro_output_basepath 'doQc/gp15/gp15_obs.mat'],         'gp15_obs');
load([sim_output_basepath 'calcMld/gp15_mld.mat'],           'gp15_mld');
load([pro_output_basepath 'doQc/ppz/gp15_ppz.mat'],          'gp15_ppz');
load([pro_output_basepath 'doQc/stations/gp15_stations.mat'],'gp15_stations', 'NUMSTAT');

%%  loop over stations
for iStat = 1 : 1 : NUMSTAT

    sn     = gp15_stations.stationNo(iStat);
    ez     = gp15_ppz.ppzDepth(gp15_ppz.stationNo == sn);
    mldVal = gp15_mld.('MLD JAK')(iStat);
    eqpVal = gp15_eqp.eqp(gp15_eqp.stationNo == sn);

    dataSn = gp15_obs(gp15_obs.stationNo == sn & gp15_obs.depth <= BTM_DEPTH, :);

    figure;
    hold('on');
    plot(dataSn.th234, dataSn.depth, '-ok', ...
         'lineWidth', 1, 'markerEdgeColor', 'k', 'markerFaceColor', 'white');
    plot(dataSn.u238, dataSn.depth, '-sk', ...
         'markerEdgeColor', 'k', 'markerFaceColor', 'w', ...
         'lineWidth', 1, 'markerSize', 5);
    yline(eqpVal,  '-k',  'lineWidth', 2);
    yline(mldVal,  '--k', 'lineWidth', 1);
    yline(ez,      '-.k', 'lineWidth', 1);
    hold('off');
    box('on');
    xlabel('\textbf{Activity  (dpm L$^{-1}$)}', 'interpreter', 'latex', 'fontSize', 20);
    ylabel('\textbf{Depth (m)}', 'interpreter', 'latex', 'fontSize', 20);
    ylim([0, 200]);
    xlim([1, 3]);
    set(gca, 'yDir', 'reverse', 'tickLabelInterpreter', 'latex', ...
        'fontSize', 16, 'fontWeight', 'bold', 'lineWidth', 1);
    set(gcf, 'units', 'inches', 'position', [2, 2, 5, 10], ...
             'paperUnits', 'inches', 'paperSize', [5, 10]);
    exportgraphics(gcf, [plot_output_basepath 'calcEqp/eqpPlots/eqpcalc' num2str(sn) '.pdf'], ...
                   'ContentType', 'vector');

end
close('all');

%%  end program
