%% regress data
%  initialize output ::
regress_coeffs = NaN(size(regions, 1), 6);  % [beta, alpha, betaVar, alphaVar, coVar, n]
deep_mean      = NaN(size(regions, 1), 3);  % [ratioSubsurface, ratioSubsurfaceStdErr, ratioSubsurfaceStd]
y_p_ppz        = cell(size(regions, 1), 1);
zstar_region   = zeros(size(regions, 1), 1);
zstar_ci       = zeros(size(regions, 1), 2);

%  set QC bounds ::
BTMDEPTH  = 400;
maxPocLsf = 3.0;
% C5: lower bound excludes LSF ratios at or near the instrumental noise floor; POC and Th-234
%     LSF measurements each carry ~5-10% relative uncertainty, so ratios below ~0.05 µmol dpm⁻¹
%     are noise-dominated and not physically meaningful. This threshold is an author-defined
%     methodological choice documented and justified in Kenyon et al. (in prep).
%     The upper bound (3.0) provides tighter domain-specific QC for the LSF fraction beyond
%     the global FLIERDATA = 30 µmol dpm⁻¹ applied in qcGp15Poc.m.
minPocLsf = 0.05;

%  compute regional mean PPZ and MLD ::
ez_region  = zeros(size(regions, 1), 1);
mld_region = zeros(size(regions, 1), 1);
for istat = 1 : size(regions, 1)
    stat_bounds      = table2array(regions(istat, 2:3));
    idx_stat         = find((gp15_stations.stationNo >= stat_bounds(1)) & (gp15_stations.stationNo <= stat_bounds(2)));
    idx_stat_mld     = find((gp15_mld.('Station No') >= stat_bounds(1)) & (gp15_mld.('Station No') <= stat_bounds(2)));
    ez_region(istat)  = mean(gp15_stations.depthPpz(idx_stat));
    mld_region(istat) = mean(gp15_mld.('MLD JAK')(idx_stat_mld));
end

%  fit and plot loop ::
figure;
tl = tiledlayout(1, size(regions, 1), 'tileSpacing', 'compact');
for istat = 1 : size(regions, 1)

    stat_bounds = table2array(regions(istat, 2:3));
    idx_stat    = find((gp15_stations.stationNo >= stat_bounds(1)) & (gp15_stations.stationNo <= stat_bounds(2)));

    %  collect all valid observations in [0, BTMDEPTH) ::
    data_all = zeros(0, 4);   % [depth, ratio, sigma, stationNo]
    for idat = 1 : length(idx_stat)
        stat    = idx_stat(idat);
        idx_obs = find((gp15_obs.stationNo == stat)                       & ...
                       (gp15_obs.depth >= 0)                              & ...
                       (gp15_obs.depth < BTMDEPTH)                        & ...
                       (~isnan(gp15_obs.pocLarge))                        & ...
                       (gp15_obs.pocTh234RatioLarge <= maxPocLsf)         & ...
                       (gp15_obs.pocTh234RatioLarge >= minPocLsf));
        if isempty(idx_obs), continue; end
        data_all = [data_all; ...
            gp15_obs.depth(idx_obs), ...
            gp15_obs.pocTh234RatioLarge(idx_obs), ...
            gp15_obs.uncertPocTh234RatioLarge(idx_obs), ...
            repmat(double(stat), numel(idx_obs), 1)];
    end
    data_all = rmmissing(data_all);

    %  fit piecewise model via profile likelihood ::
    [fittedAlpha, fittedBeta, fittedZStar, paramCov, zStarCI, ~, nBinMeans] = ...
        fitPiecewiseRatio(data_all(:,1), data_all(:,2), data_all(:,3), ...
                          data_all(:,4), mld_region(istat));

    %  store regression coefficients ::
    regress_coeffs(istat, :) = [fittedBeta, fittedAlpha, ...
                                 paramCov(2,2), paramCov(1,1), paramCov(1,2), nBinMeans];
    zstar_region(istat)  = fittedZStar;
    zstar_ci(istat, :)   = zStarCI;

    %  deep-zone constant: value at zStar by continuity ::
    ratioAtZStar = fittedAlpha + fittedBeta * fittedZStar;
    xiZStar      = [1, fittedZStar];
    sigAtZStar   = sqrt(max(xiZStar * paramCov * xiZStar', 0));
    deep_mean(istat, :) = [ratioAtZStar, sigAtZStar, sigAtZStar];

    %  evaluate piecewise fit on plotting grid ::
    depth_eval = linspace(0, BTMDEPTH, 1000)';
    zActive    = max(min(depth_eval, fittedZStar), mld_region(istat));
    y_fit_eval = fittedAlpha + fittedBeta * zActive;
    sigFit     = sqrt(max(paramCov(1,1) + zActive.^2 .* paramCov(2,2) ...
                          + 2 .* zActive .* paramCov(1,2), 0));
    y_p_ppz{istat} = [y_fit_eval + sigFit, y_fit_eval - sigFit];  % 1000 x 2

    %  plot ::
    nexttile();
    hold('on');
    hCi  = fill([y_p_ppz{istat}(:,1); flipud(y_p_ppz{istat}(:,2))], ...
                [depth_eval;           flipud(depth_eval)], ...
                [0.7 0.7 0.7], 'faceAlpha', 0.4, 'edgeColor', 'none');
    hDat = errorbar(data_all(:,2), data_all(:,1), data_all(:,3), data_all(:,3), 'horizontal', ...
                    'marker', 'o', 'markerSize', 10, ...
                    'markerEdgeColor', 'k', 'markerFaceColor', 'w', ...
                    'color', 'k', 'lineStyle', 'none', 'lineWidth', 1);
    hDat.CapSize = 0;
    hFit = plot(y_fit_eval, depth_eval, 'k', 'lineWidth', 2);
    hMld = yline(mld_region(istat), 'k--', 'lineWidth', 1);
    hZs  = yline(fittedZStar,       'r:',  'lineWidth', 1.5);
    hEz  = yline(ez_region(istat),  'k-.', 'lineWidth', 1);
    h100 = yline(100,               'k:',  'lineWidth', 1);
    title(region_names(istat), 'interpreter', 'latex');
    set(gca, 'box', 'on', 'yDir', 'reverse', ...
        'xLim', [0, 3], 'yLim', [0, BTMDEPTH], 'tickLabelInterpreter', 'latex');
    if istat == size(regions, 1)
        legend([hDat, hFit, hCi, hMld, hZs, hEz, h100], ...
               {'Data $\pm\,\sigma$', 'Fit', '$\pm 1\sigma$', 'MLD', '$z^*$', 'PPZ', '100 m'}, ...
               'interpreter', 'latex', 'location', 'southEast');
    end

    %  accumulate plot output structs ::
    if istat == 1
        plot_output_data.data_region        = repmat(istat, size(data_all,1), 1);
        plot_output_data.depth              = data_all(:,1);
        plot_output_data.data               = data_all(:,2);
        plot_output_fit.fit_region          = repmat(istat, 1000, 1);
        plot_output_fit.fit_depth           = depth_eval;
        plot_output_fit.fit_data            = y_fit_eval;
        plot_output_fit.fit_data_error_low  = y_p_ppz{istat}(:,1);
        plot_output_fit.fit_data_error_high = y_p_ppz{istat}(:,2);
    else
        plot_output_data.data_region        = [plot_output_data.data_region;        repmat(istat, size(data_all,1), 1)];
        plot_output_data.depth              = [plot_output_data.depth;              data_all(:,1)];
        plot_output_data.data               = [plot_output_data.data;               data_all(:,2)];
        plot_output_fit.fit_region          = [plot_output_fit.fit_region;          repmat(istat, 1000, 1)];
        plot_output_fit.fit_depth           = [plot_output_fit.fit_depth;           depth_eval];
        plot_output_fit.fit_data            = [plot_output_fit.fit_data;            y_fit_eval];
        plot_output_fit.fit_data_error_low  = [plot_output_fit.fit_data_error_low;  y_p_ppz{istat}(:,1)];
        plot_output_fit.fit_data_error_high = [plot_output_fit.fit_data_error_high; y_p_ppz{istat}(:,2)];
    end

end

%  set labels ::
xlabel(tl, 'LSF POC:$^{234}$Th [$\mu$mol dpm$^{-1}$]', 'interpreter', 'latex');
ylabel(tl, 'Depth [m]', 'interpreter', 'latex');
set(gcf, 'position', [0, 0, 1200, 400]);
exportgraphics(gcf, [plot_output_basepath 'regressRegionalRatio/regressionPocLsf.pdf'], 'ContentType', 'vector');

%% save final result
%  make table ::
ratioPocLsf = array2table([regress_coeffs, deep_mean, mld_region, zstar_region, zstar_ci, ez_region]);
ratioPocLsf = [region_names', ratioPocLsf];
ratioPocLsf.Properties.VariableNames = {'region', 'ratioBeta', 'ratioAlpha', 'betaVar', 'alphaVar', 'coVar', 'n', ...
                                          'ratioSubsurface', 'ratioSubsurfaceStdErr', 'ratioSubsurfaceStd', ...
                                          'mldMean', 'zstar', 'zstarCI_lo', 'zstarCI_hi', 'ppzmean'};
plot_output_data = struct2table(plot_output_data);
plot_output_fit  = struct2table(plot_output_fit);

%  save ::
save([sim_output_basepath 'regressRegionalRatio/ratioPocLsf.mat'], 'ratioPocLsf');
save([sim_output_basepath 'regressRegionalRatio/plotting/plot_output_data.mat'], 'plot_output_data');
save([sim_output_basepath 'regressRegionalRatio/plotting/plot_output_fit.mat'],  'plot_output_fit');
writetable(ratioPocLsf,       [sim_output_basepath 'regressRegionalRatio/ratioPocLsf.xlsx']);
writetable(plot_output_data,  [sim_output_basepath 'regressRegionalRatio/plotting/plot_output_data.xlsx']);
writetable(plot_output_fit,   [sim_output_basepath 'regressRegionalRatio/plotting/plot_output_fit.xlsx']);

%% end subroutine
