%%  sensitivity: spline smoothing parameter for vertical Th-234 gradient
%   refits the gradient spline for pTest = {0.45, 0.9, 1.0} at one
%   representative equatorial and one representative subtropical station.
%   runs after doVertGradCalc; inside calcVertGrad.m. ::

%%  define test values
pTest      = [0.45, 0.9, 1.0];
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

    iStat   = testStatIdxs(iPlot);
    sn      = testStatNos(iPlot);
    ez      = gp15_ppz.ppzDepth(gp15_ppz.stationNo == sn);
    mldVal  = gp15_mld.('MLD JAK')(iStat);
    rootsFt = gp15_eqp.eqp(gp15_eqp.stationNo == sn);
    zBot    = min(rootsFt, ez);
    statLat = gp15_stations.latitude(iStat);

    %   extract profile ::
    x   = gp15_obs.depth(gp15_obs.stationNo == sn);
    yth = gp15_obs.th234(gp15_obs.stationNo == sn);
    nth = gp15_obs.uncertTh234(gp15_obs.stationNo == sn);
    n   = length(x);

    %   raw finite differences (same for all p values) ::
    gradTh      = zeros(n, 1);
    deltaGradTh = zeros(n, 1);
    for iDepth = 1 : 1 : n
        if iDepth == 1
            th1 = yth(1);   th2 = yth(2);
            dt1 = nth(1);   dt2 = nth(2);
            z1  = x(1);     z2  = x(2);
        elseif iDepth == n
            th1 = yth(n-1); th2 = yth(n);
            dt1 = nth(n-1); dt2 = nth(n);
            z1  = x(n-1);   z2  = x(n);
        else
            th1 = yth(iDepth-1); th2 = yth(iDepth+1);
            dt1 = nth(iDepth-1); dt2 = nth(iDepth+1);
            z1  = x(iDepth-1);   z2  = x(iDepth+1);
        end
        B                   = abs(z2 - z1);
        gradTh(iDepth)      = (th2 - th1) / (z2 - z1);
        deltaGradTh(iDepth) = (dt1 / B) ^ 2 + (dt2 / B) ^ 2;   % variance
    end

    %   spline augmentation: boundary anchors at mldVal and zBot (same as doVertGradCalc) ::
    weights  = 1 ./ deltaGradTh;
    W_anchor = max(weights);
    x_aug    = [x;        mldVal;   zBot    ];
    gTh_aug  = [gradTh;   0;        0       ];
    wts_aug  = [weights;  W_anchor; W_anchor];
    [x_aug, sortIdx] = sort(x_aug);
    gTh_aug  = gTh_aug(sortIdx);
    wts_aug  = wts_aug(sortIdx);

    %   spline fit for each p value; apply same active-zone mask as doVertGradCalc ::
    gradHats = NaN(n, nTest);
    for iP = 1 : 1 : nTest
        gh = csaps(x_aug, gTh_aug, pTest(iP), x, wts_aug);
        gh(x < mldVal | x > zBot) = 0;
        gradHats(:, iP) = gh;
    end

    %   report gradient at PPZ for each p value (proxy for 2D flux correction sensitivity) ::
    [~, iPpz] = min(abs(x - ez));
    fprintf('\n%s (stn %d, lat %.1f deg N): gradient at PPZ (%.0f m)\n', ...
            testStatLabels{iPlot}, sn, statLat, ez);
    for iP = 1 : 1 : nTest
        fprintf('  p = %.2f:  grad = %+.4f dpm L-1 m-1\n', pTest(iP), gradHats(iPpz, iP));
    end

    %   plot ::
    figure;
    hold('on');
    plot(gradTh, x, 'o', ...
         'color', [0.6, 0.6, 0.6], ...
         'markerSize', 10, ...
	 'markerFaceColor', 'w', ...
         'markerEdgeColor', 'k', ...
         'displayName', 'Finite differences');
    for iP = 1 : 1 : nTest
        plot(gradHats(:, iP), x, lineStyles{iP}, ...
             'color', 'k', ...
             'lineWidth', 2, ...
             'displayName', ['$p = ' num2str(pTest(iP)) '$']);
    end
    xline(0, '-', 'color', [0.5, 0.5, 0.5], 'lineWidth', 0.5, 'handleVisibility', 'off');
    yline(mldVal, '--', 'color', [0.00, 0.45, 0.74], 'lineWidth', 1, ...
          'label', 'MLD', 'labelHorizontalAlignment', 'left', ...
          'handleVisibility', 'off');
    yline(zBot, '--', 'color', [0.85, 0.33, 0.10], 'lineWidth', 1, ...
          'label', '$z_\mathrm{bot}$', 'labelHorizontalAlignment', 'left', ...
          'interpreter', 'latex', 'handleVisibility', 'off');
    hold('off');
    box('on');
    xlabel('\textbf{Gradient (dpm L$^{-1}$ m$^{-1}$)}', 'interpreter', 'latex', 'fontSize', 20);
    ylabel('\textbf{Depth (m)}', 'interpreter', 'latex', 'fontSize', 20);
    legend('location', 'southEast', 'interpreter', 'latex', 'fontSize', 14);
    title(['\textbf{Station ' num2str(sn) ' (' sprintf('%.1f', statLat) '$^{\circ}$N)}'], ...
          'interpreter', 'latex', 'fontSize', 18);
    ylim([0, 400]);
    set(gca, 'yDir', 'reverse', 'tickLabelInterpreter', 'latex', ...
        'fontSize', 16, 'fontWeight', 'bold', 'lineWidth', 1);
    set(gcf, 'units', 'inches', 'position', [2, 2, 5, 10], ...
             'paperUnits', 'inches', 'paperSize', [5, 10]);
    exportgraphics(gcf, ...
                   [plot_output_basepath 'calcVertGrad/sensitivity_splineGrad_stn' ...
                    num2str(sn) '.pdf'], 'ContentType', 'vector');

end

%%  end subroutine
