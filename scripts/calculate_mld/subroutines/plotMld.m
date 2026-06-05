%% plotMld - diagnostic plots of mixed layer depth profiles
%  standalone; run after calcMld.m has saved gp15_mld.mat.
%  recomputes the density spline and mixed-layer statistics from gp15_ctd
%  and the subjective MLD guesses; MLD values come from the saved output.
%--------------------------------------------------------------------------
setupModel;
setModelCoefficients;
configureCalcMld;

%%  load saved outputs and raw inputs
load([sim_output_basepath 'calcMld/gp15_mld.mat'], 'gp15_mld');
load([pro_output_basepath 'doQc/gp15/gp15_ctd.mat'], 'gp15_ctd');
load([pro_output_basepath 'doQc/stations/gp15_stations.mat'], 'gp15_stations', 'NUMSTAT');

%   subjective MLD guesses (needed to define mixed-layer density band) ::
mldCalcs = readtable([input_basepath 'mld/mld_calcs_jak.xlsx'], 'variableNamingRule', 'preserve');

%%  loop over stations
for iStat = 1 : 1 : NUMSTAT

    sn       = gp15_stations.stationNo(iStat);
    statCast = gp15_stations.castNo(iStat);

    %   MLD values from saved output ::
    mld_JAK_stat = gp15_mld.('MLD JAK')(iStat);
    mld_dBM_stat = gp15_mld.('MLD dBM')(iStat);
    mld_BW_stat  = [gp15_mld.('MLD BW 01')(iStat), ...
                    gp15_mld.('MLD BW 05')(iStat), ...
                    gp15_mld.('MLD BW 125')(iStat)];

    %   density profile for this station/cast ::
    iCtdStatCast  = find(gp15_ctd.stationNo == sn & gp15_ctd.castNo == statCast);
    ctdStatCast   = gp15_ctd(iCtdStatCast, :);
    mldGuess      = mldCalcs.mld_JAK(iStat);
    iDepthPlot    = find(ctdStatCast.ctd_depth >= mldGuess + ADD_DEPTH);
    if isempty(iDepthPlot)
        endDepthPlot = height(ctdStatCast);
    else
        endDepthPlot = iDepthPlot(1);
    end
    ctdPlot       = ctdStatCast(1 : endDepthPlot, :);
    potDensProf   = ctdPlot.ctd_potentialDensity;
    potDensDepthProf = ctdPlot.ctd_depth;

    %   density spline for visualization ::
    potDensSpline    = csaps(potDensDepthProf, potDensProf, P_SPLINE_DENSITY, [], ones(size(potDensProf)));
    potDensPlotDepth = 0 : 0.1 : max(potDensDepthProf);
    potDensPlot      = ppval(potDensSpline, potDensPlotDepth);

    %   mixed-layer density band (needed for xline markers) ::
    potDensMld        = potDensProf(potDensDepthProf <= mldGuess);
    meanMld_stat      = mean(potDensMld, 'all');
    stdMld_stat       = std(potDensMld, 0, 'all');

    %   plot ::
    figure;

    hold('on');
    plot(potDensPlot, potDensPlotDepth, '--k', 'lineWidth', 1);
    scatter(potDensProf, potDensDepthProf, 50, 's', ...
            'lineWidth', 1, 'markerEdgeColor', 'k', 'markerFaceColor', 'white');
    yline(mld_JAK_stat, '-k', 'lineWidth', 2);
    yline(mld_dBM_stat, '-k', 'lineWidth', 1);
    yline(mld_BW_stat(1), '--k', 'lineWidth', 1);
    yline(mld_BW_stat(2), '-.k', 'lineWidth', 1);
    yline(mld_BW_stat(3), ':k', 'lineWidth', 1);
    xline(meanMld_stat + (NUM_STD * stdMld_stat), ':k', 'lineWidth', 0.5);
    xline(meanMld_stat - (NUM_STD * stdMld_stat), ':k', 'lineWidth', 0.5);
    hold('off');

    legend('Potential density $\rho$ smoothed spline', ...
	   'CTD potential density', ...
	   'Subject-objective $\langle\rho(z)\rangle_\textrm{MLD}$ (Pickart et al.)', ...
	   'Fixed $\Delta T = 0.2$ (De Boyer Montegut, 2004)', ...
	   'Fixed $\Delta \rho = 0.01$ (Bishop and Wood, 2009)', ...
	   'Fixed $\Delta \rho = 0.05$ (Bishop and Wood, 2009)', ...
	   'Fixed $\Delta \rho = 0.125$ (Bishop and Wood, 2009)', ...
	   ['$\langle\rho(z)\rangle_\textrm{MLD} \:\pm\: ' num2str(NUM_STD), '\sigma$'], ...
	   'fontSize', 16, 'fontWeight', 'bold', 'interpreter', 'latex', 'box', 'off', 'location', 'eastoutside');

    box('on');
    xlabel('\textbf{Potential Density  (kg m$^{-3}$)}', 'interpreter', 'latex', 'fontSize', 20);
    ylabel('\textbf{Depth (m)}', 'interpreter', 'latex', 'fontSize', 20);
    ylim([0, 80]);
    set(gca, 'yDir', 'reverse', 'tickLabelInterpreter', 'latex', ...
        'fontSize', 16, 'fontWeight', 'bold', 'lineWidth', 1);
    set(gcf, 'units', 'inches', 'position', [2, 2, 12.5, 10], ...
             'paperUnits', 'inches', 'paperSize', [12.5, 10]);
    exportgraphics(gcf, [plot_output_basepath 'calcMld/mldPlots/mldcalc' num2str(sn) '.pdf'], ...
                   'ContentType', 'vector');

end
close('all');

%%  end program
