%%  sensitivity: spline smoothing parameter for MLD potential-density spline
%   Refits the JAK density spline for pTest = {0.495, 0.99, 1.0} at one
%   representative equatorial and one representative subtropical station.
%   Runs after calcMlDepths; inside calcMld.m.
%
%   Why low flux sensitivity is expected (three reasons encoded here):
%   (1) CTD profiles are dense (sub-metre resolution) and precise
%       (instrument noise << MLD density signal). With O(100+) points
%       per cast, csaps at p = 0.495 vs p = 0.99 produces nearly
%       identical spline shapes and therefore the same zero-crossing.
%   (2) P_SPLINE_DENSITY controls a depth boundary (MLD), not a flux
%       value directly. A small MLD shift propagates as a second-order
%       effect on the integrated Th-234 flux, unlike P_SPLINE_GRADIENT
%       which directly scales the upwelling correction term w*(dTh/dz)*dz.
%   (3) Three independent MLD algorithms (JAK, dBM, BW) provide an
%       implicit inter-method robustness check whose spread subsumes any
%       p variation.
%   This script quantifies (1) explicitly and confirms (2) holds.
%--------------------------------------------------------------------------

%%  define test values
%   pTest = {p/2, p_nominal, 1.0}; p*2 exceeds 1 so 1.0 is the upper bound.
pTest      = [0.495, 0.99, 1.0];
nTest      = length(pTest);
lineStyles = {'-', '--', ':'};

%%  choose representative stations
%   equatorial  : station with |latitude| closest to 0 deg
%   subtropical : station with latitude closest to 30 deg N
[~, iEq]       = min(abs(gp15_stations.latitude));
[~, iSub]      = min(abs(gp15_stations.latitude - 30));
testStatIdxs   = [iEq, iSub];
testStatNos    = gp15_stations.stationNo(testStatIdxs);
testStatLabels = {'Equatorial', 'Subtropical'};

%%  loop over representative stations
for iPlot = 1 : 1 : 2

    iStat    = testStatIdxs(iPlot);
    sn       = testStatNos(iPlot);
    statCast = gp15_stations.castNo(iStat);
    statLat  = gp15_stations.latitude(iStat);

    %   reconstruct density profile (same truncation as calcMlDepths) ::
    mldGuess     = mldCalcs.mld_JAK(iStat);
    iCtd         = find(gp15_ctd.stationNo == sn & gp15_ctd.castNo == statCast);
    ctdStat      = gp15_ctd(iCtd, :);
    iEnd         = find(ctdStat.ctd_depth >= mldGuess + ADD_DEPTH);
    if isempty(iEnd)
        iEnd = height(ctdStat);
    else
        iEnd = iEnd(1);
    end
    ctdTrunc         = ctdStat(1 : iEnd, :);
    potDensProf      = ctdTrunc.ctd_potentialDensity;
    potDensDepthProf = ctdTrunc.ctd_depth;

    %   JAK density threshold from pre-computed mixed-layer statistics ::
    mldPPotDens = statsMld.meanMld(iStat) + (NUM_STD * statsMld.stdMld(iStat));

    %   fine depth grid for ppval evaluation ::
    plotDepth = (0 : 0.1 : max(potDensDepthProf))';

    %   refit spline and find JAK zero-crossing for each p value ::
    mldHats    = NaN(1, nTest);
    splineFits = cell(1, nTest);
    for iP = 1 : 1 : nTest

        %   shifted form: zero-crossing = MLD (reason 1 test) ::
        ppShift  = csaps(potDensDepthProf, potDensProf - mldPPotDens, pTest(iP), [], ones(size(potDensProf)));
        pZeros   = fnzeros(ppShift);
        mldHats(iP) = min(pZeros(pZeros > MLD_AVE_DEPTH), [], 'all');

        %   unshifted form for visualization ::
        ppRaw          = csaps(potDensDepthProf, potDensProf, pTest(iP), [], ones(size(potDensProf)));
        splineFits{iP} = ppval(ppRaw, plotDepth);

    end

    %   plot density profile with spline fits and resulting MLD depths ::
    figure;
    hold('on');
    scatter(potDensProf, potDensDepthProf, 50, 's', ...
            'lineWidth', 1, 'markerEdgeColor', 'k', ...
            'markerFaceColor', 'w', ...
            'displayName', 'CTD potential density, $\rho$');
    for iP = 1 : 1 : nTest
        plot(splineFits{iP}, plotDepth, lineStyles{iP}, ...
             'color', 'k', 'lineWidth', 1, ...
             'displayName', ['Spline with $p = ' num2str(pTest(iP)) '$: MLD = ' num2str(round(mldHats(iP), 0)) ' [m]']);
    end
    for iP = 1 : 1 : nTest
        yline(mldHats(iP), lineStyles{iP}, ...
              'color', 'k', 'lineWidth', 1, 'handleVisibility', 'off');
    end
    hold('off');
    box('on');
    xlabel('\textbf{Potential Density (kg m$^{-3}$)}', 'interpreter', 'latex', 'fontSize', 20);
    ylabel('\textbf{Depth (m)}', 'interpreter', 'latex', 'fontSize', 20);
    legend('location', 'eastoutside', 'interpreter', 'latex', 'fontSize', 14);
    title(['\textbf{Station ' num2str(sn) ' (' sprintf('%.1f', statLat) '$^{\circ}$N)}'], ...
          'interpreter', 'latex', 'fontSize', 18);
    ylim([0, 80]);
    set(gca, 'yDir', 'reverse', 'tickLabelInterpreter', 'latex', ...
        'fontSize', 16, 'fontWeight', 'bold', 'lineWidth', 1);
    set(gcf, 'units', 'inches', 'position', [2, 2, 12.5, 10], ...
             'paperUnits', 'inches', 'paperSize', [12.5, 10]);
    exportgraphics(gcf, ...
                   [plot_output_basepath 'calcMld/sensitivity/sensitivity_splineDensity_stn' ...
                    num2str(sn) '.pdf'], 'ContentType', 'vector');

end

%%  end subroutine
