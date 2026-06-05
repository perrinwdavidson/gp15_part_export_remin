%% load data
%   all products (ecco and mercator) are loaded from their calcUpwell netCDF
%   output, which contains wSpatAve — the spatially and 35-d temporally
%   averaged field. ecco previously loaded raw ecco_wo.mat here; it now
%   enters via the same averaged pipeline as all other products so that
%   kriging training data is consistent in magnitude and smoothing. ::
fileName = [sim_output_basepath 'calcUpwell/w_' dataProduct '.nc'];
X0 = ncread(fileName, 'longitude');
Y0 = ncread(fileName, 'latitude');
Z0 = ncread(fileName, 'depth');
T0 = ncread(fileName, 'time');
V0 = ncread(fileName, 'w');

%%  data prep
%%% make coordinates ::
[lon, lat, depth, time] = ndgrid(X0, Y0, Z0, T0);

%%% get data ::
X = lon((X0 >= pacificLon(1)) & (X0 <= pacificLon(2)), ...
        (Y0 >= pacificLat(1)) & (Y0 <= pacificLat(2)), ...
        (Z0 >= pacificDepth(1)) & (Z0 <= pacificDepth(2)), ...
        :);
Y = lat((X0 >= pacificLon(1)) & (X0 <= pacificLon(2)), ...
        (Y0 >= pacificLat(1)) & (Y0 <= pacificLat(2)), ...
        (Z0 >= pacificDepth(1)) & (Z0 <= pacificDepth(2)), ...
        :);
Z = depth((X0 >= pacificLon(1)) & (X0 <= pacificLon(2)), ...
          (Y0 >= pacificLat(1)) & (Y0 <= pacificLat(2)), ...
          (Z0 >= pacificDepth(1)) & (Z0 <= pacificDepth(2)), ...
          :);
T = time((X0 >= pacificLon(1)) & (X0 <= pacificLon(2)), ...
         (Y0 >= pacificLat(1)) & (Y0 <= pacificLat(2)), ...
         (Z0 >= pacificDepth(1)) & (Z0 <= pacificDepth(2)), ...
         :);
V = V0((X0 >= pacificLon(1)) & (X0 <= pacificLon(2)), ...
       (Y0 >= pacificLat(1)) & (Y0 <= pacificLat(2)), ...
       (Z0 >= pacificDepth(1)) & (Z0 <= pacificDepth(2)), ...
       :);

%%  set initial GP hyperparameters
%   physically motivated priors; re-optimised by ordinary_kriging (FitHp=true, NRestarts=3).
%   the ecco_pacific.nc save block and calculate_decorrelation.jl pre-step that previously
%   derived these from data have been removed — fitting is done internally by the kriging. ::
lxInit = 500e3;  % [m] horizontal scale (~500 km, order-of-magnitude Pacific basin prior)
lt     = 30.0;   % [d] temporal scale (~1 month, intraseasonal upwelling prior)

%% previously: detrend ecco and save for external decorrelation analysis
% if strcmp(dataProduct, 'ecco')
%     Vtilde = V;
%     for i = 1 : 1 : size(V, 4)
%         Vtilde(:, :, :, i) = V(:, :, :, i) - nanmean(V(:, :, :, i), 'all');
%     end
%     pacificData = rmmissing([X(:), Y(:), T(:), V(:), Vtilde(:)]);
%     fname = [sim_output_basepath 'interpData/decorrelation/ecco_pacific.nc'];
%     delete(fname);
%     vname = 'pacific_wo';
%     nccreate(fname, vname, 'dimensions', {'observations', size(pacificData, 1), 'variables', size(pacificData, 2)}, 'fillValue', 'disable');
%     ncwrite(fname, vname, pacificData);
% end

%% previously: load decorrelation scale CSVs produced by calculate_decorrelation.jl
% lxRaw = readmatrix([sim_output_basepath 'interpData/decorrelation/spatial_decorrelation.csv']);
% ltRaw = readmatrix([sim_output_basepath 'interpData/decorrelation/temporal_decorrelation.csv']);
% lxInit = lxRaw(2, 1) * 1000;  % [km] -> [m]
% lt     = round(ltRaw(2, 1),  1 + int64(floor(log10(ltRaw(2, 2)))), 'significant');  % [d]

%%  build flat sample array (all non-missing, non-zero model grid points)
%   columns: [lon [deg, 0-360], lat [deg], depth [m], time [d], w [m s-1]] ::
Xs = rmmissing([X(:), Y(:), Z(:), T(:), V(:)]);
Xs(Xs(:, 5) == 0, :) = [];  % remove zero-velocity (land / unfilled) points

%%  fit depth-binned temporal preprocessor on training w
%   fits annual + semi-annual harmonic + linear trend per depth bin;
%   returns anomaly observations and path to saved preprocessor for inversion ::
[anomV, preproc_path] = fit_temporal_preprocessor(Xs(:, 2), Xs(:, 1), Xs(:, 3), Xs(:, 5), Xs(:, 4), ...
                                                   'DepthBinEdges', [0, 100, 200, 400, Inf], ...
                                                   'RandomState',   7);

%%  build query arrays — one row per (station, observation depth) pair
query_lat   = [];
query_lon   = [];
query_depth = [];
query_time  = [];
query_sn    = [];

for i = 1 : 1 : NUMSTAT

	sn = gp15_stations.stationNo(i);
	zq = gp15_obs.depth((gp15_obs.stationNo == sn) & (gp15_obs.depth <= BTMDEPTH));
	tq = datenum(gp15_stations.date(i));
	nz = length(zq);

	query_lat   = [query_lat;   repmat(gp15_stations.latitude(i),        nz, 1)];
	query_lon   = [query_lon;   repmat(gp15_stations.longitude(i) + 360, nz, 1)];
	query_depth = [query_depth; zq];
	query_time  = [query_time;  repmat(tq,                               nz, 1)];
	query_sn    = [query_sn;    repmat(sn,                               nz, 1)];

end

%%  spatiotemporal ordinary kriging — single call over all (station, depth) query points
[wAnom, wStd] = ordinary_kriging(Xs(:, 2), Xs(:, 1), Xs(:, 3), anomV, ...
                                  query_lat, query_lon, query_depth, ...
                                  'TrainTime',          Xs(:, 4), ...
                                  'QueryTime',          query_time, ...
                                  'GpHorizontalScaleM', lxInit, ...
                                  'GpVerticalScaleM',   100, ...
                                  'GpTimeScale',        lt, ...
                                  'GpTimeScaleMin',     1.0, ...
                                  'GpTimeScaleMax',     365.25, ...
                                  'MaxK',               500, ...
                                  'MinK',               10, ...
                                  'FitHp',              true, ...
				  'Variance', 		1.0, ...
				  'NoiseVariance', 	0.1, ...
                                  'Kernel',             'matern32', ...
                                  'NRestarts',          3, ...
                                  'RandomState',        7);

%%  invert temporal preprocessing to recover physical units [m s-1]
w_est = inverse_temporal_preprocess(wAnom, query_time, query_depth, preproc_path);

%%  package output table
%   wErr = wStd (ordinary-kriging standard deviation); wVar = wStd^2.
%   wN = 1 for all rows (kriging gives a single BLUP per query point, not an average).
%   column order matches the prior format for backward compatibility with downstream loaders. ::
gp15_w = table(query_sn, ...
               query_lon - 360, ...
               query_lat, ...
               query_depth, ...
               query_time, ...
               w_est, ...
               wStd .^ 2, ...
               wStd, ...
               ones(length(query_sn), 1), ...
               'VariableNames', {'stationNo', 'longitude', 'latitude', 'depth', 'time', 'w', 'wVar', 'wErr', 'wN'});

%%  diagnostic plots
%%% extract slices at 100 m and PPZ ::
w100 = gp15_w(gp15_w.depth == 100, :);
[~, idx100, ~] = unique(w100.latitude);
w100 = w100(idx100, :);

% diagnostic: old: exact float equality — silently returns empty wPpz if kriging depths have any rounding
% [~, idxPpz0] = intersect(gp15_w.depth, gp15_stations.depthPpz);
% wPpz = gp15_w(idxPpz0, :);
idxPpz0 = false(height(gp15_w), 1);
for iPpz = 1 : 1 : NUMSTAT
	snPpz  = gp15_stations.stationNo(iPpz);
	ppzDep = gp15_stations.depthPpz(iPpz);
	idxPpz0 = idxPpz0 | ((gp15_w.stationNo == snPpz) & (abs(gp15_w.depth - ppzDep) <= 1e-6));
end
wPpz = gp15_w(idxPpz0, :);
[~, idxPpz, ~] = unique(wPpz.latitude);
wPpz = wPpz(idxPpz, :);

%%% scatter: latitude × depth, coloured by w ::
figure;
tl = tiledlayout(2, 1, 'tileSpacing', 'compact');

nexttile();
scatter(gp15_w.latitude, gp15_w.depth, 200, gp15_w.w * SEC2DAY, 'filled');
xlabel('\textbf{Latitude [deg N.]}', 'interpreter', 'latex');
ylabel('\textbf{Depth [m]}', 'interpreter', 'latex');
set(gca, 'yDir', 'reverse', 'tickLabelInterpreter', 'latex', 'fontSize', 16, 'fontWeight', 'bold', 'lineWidth', 1);
title(['\textbf{PMT Upwelling Velocity (' upper(dataProduct) ')}'], 'interpreter', 'latex', 'fontSize', 20);
cb = colorbar;
ylabel(cb, '$w$ [m d$^{-1}$]', 'interpreter', 'latex');
set(cb, 'tickLabelInterpreter', 'latex', 'fontSize', 16, 'fontWeight', 'bold', 'lineWidth', 1);
caxis([-5 5]);

nexttile();
scatter(gp15_w.latitude, gp15_w.depth, 200, gp15_w.wErr * SEC2DAY, 'filled');
xlabel('\textbf{Latitude [deg N.]}', 'interpreter', 'latex');
ylabel('\textbf{Depth [m]}', 'interpreter', 'latex');
set(gca, 'yDir', 'reverse', 'tickLabelInterpreter', 'latex', 'fontSize', 16, 'fontWeight', 'bold', 'lineWidth', 1);
title(['\textbf{PMT Upwelling Velocity Standard Deviation (' upper(dataProduct) ')}'], 'interpreter', 'latex', 'fontSize', 20);
cb = colorbar;
ylabel(cb, '$\sigma_{\hat{w}}$ [m d$^{-1}$]', 'interpreter', 'latex');
set(cb, 'tickLabelInterpreter', 'latex', 'fontSize', 16, 'fontWeight', 'bold', 'lineWidth', 1);
caxis([-0.5 0.5]);

set(gcf, 'position', [0, 0, 1000, 1000]);
exportgraphics(gcf, [plot_output_basepath 'interpData/upwelling/' dataProduct '_transect.pdf'], 'ContentType', 'vector');

%%% line plot: w at 100 m and PPZ vs latitude ::
figure;
hold('on');

plot(w100.latitude, w100.w * SEC2DAY, '--k', 'lineWidth', 1.5)
err1 = errorbar(w100.latitude, w100.w * SEC2DAY, w100.wErr * SEC2DAY, w100.wErr * SEC2DAY, '--k', 'lineWidth', 1.5);
err1.Color = [0 0 0];
err1.LineStyle = 'none';

plot(wPpz.latitude, wPpz.w * SEC2DAY, '-k', 'lineWidth', 1.5)
err2 = errorbar(wPpz.latitude, wPpz.w * SEC2DAY, wPpz.wErr * SEC2DAY, wPpz.wErr * SEC2DAY, '-k', 'lineWidth', 1.5);
err2.Color = [0 0 0];
err2.LineStyle = 'none';

yline(0, ':k', 'lineWidth', 1.5);
hold('off');

xlabel('\textbf{Latitude (deg. N)}', 'interpreter', 'latex', 'fontSize', 20);
ylabel('\textbf{Upwelling velocity, $w$ (m day$^{-1}$)}', 'interpreter', 'latex', 'fontSize', 20);
legend('$w_{100}$', '', '$w_{PPZ}$', '$\epsilon_{w}$', 'interpreter', 'latex', 'fontSize', 16, 'location', 'eastOutside');

set(gca, 'box', 'on', 'tickLabelInterpreter', 'latex', 'fontSize', 16, 'fontWeight', 'bold', 'lineWidth', 1);
set(gcf, 'units', 'inches', 'position', [1.5, 1.5, 16, 7], 'paperUnits', 'inches', 'paperSize', [16, 7]);
exportgraphics(gcf, [plot_output_basepath 'interpData/upwelling/' dataProduct '_transect_depths.pdf'], 'ContentType', 'vector');

%%  save data
%   ecco: variable gp15_w  in w_ecco.mat
%   mercator: variable gp15_wAve in w_<product>Ave.mat
%   (naming preserved for backward compatibility with compareModels and loadCollateData) ::
if strcmp(dataProduct, 'ecco')

	writetable(gp15_w, [sim_output_basepath 'interpData/upwelling/w_' dataProduct '.xlsx']);
	save([sim_output_basepath 'interpData/upwelling/w_' dataProduct '.mat'], 'gp15_w');

else

	gp15_wAve = gp15_w;
	writetable(gp15_wAve, [sim_output_basepath 'interpData/upwelling/w_' dataProduct 'Ave.xlsx']);
	save([sim_output_basepath 'interpData/upwelling/w_' dataProduct 'Ave.mat'], 'gp15_wAve');

end

%%  done with subroutine
