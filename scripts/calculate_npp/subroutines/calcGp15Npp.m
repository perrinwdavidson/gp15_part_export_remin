%%  loop through all measurement locations and interpolate 
%   make sample array ::
Xs = rmmissing([X(:), ...  % longitude [deg]
		Y(:), ...  % latitude [deg] 
	        T(:), ...  % time [d]
	        V(:)]);    % NPP [mmolC]

%   fit temporal pre-processor ::
[anomV, preproc] = fit_temporal_preprocessor(Xs(:, 2), ...
				             Xs(:, 1), ...
				             [], ...
			     	             Xs(:, 4), ...
				             Xs(:, 3));

%   get query points ::
sn = gp15_stations.stationNo;
xq = gp15_stations.longitude;
yq = gp15_stations.latitude;
tq = datenum(gp15_stations.date);

%   build 35-day trailing window query points for time-matched EzRatio NPP ::
%   the Th-234 inventory integrates over ~35 d prior to sampling (mean lifetime
%   = 24.101/ln2 ≈ 34.8 d); NPP must be averaged over the same trailing window
%   for a temporally consistent EzRatio numerator and denominator
nDaysWindow = 35;   % Th-234 mean lifetime rounded to nearest integer [d]
nStations   = length(sn);
dayOffsets  = (-(nDaysWindow - 1) : 0)';   % 35×1: [-34, -33, ..., 0]

%   tqWindowMat is nStations×nDaysWindow; row i is the 35 daily query times
%   for station i covering [t_sample-34, ..., t_sample]
tqWindowMat  = tq(:) + dayOffsets';   % nStations × nDaysWindow

%   flatten station-major (all 35 days of station 1, then station 2, ...) by
%   transposing to nDaysWindow×nStations then column-major reshape
tqWindowFlat = reshape(tqWindowMat.', [], 1);   % (nDaysWindow*nStations) × 1
yqWindowFlat = repelem(yq(:), nDaysWindow);      % same lat repeated nDaysWindow× per station
xqWindowFlat = repelem(xq(:), nDaysWindow);      % same lon repeated nDaysWindow× per station

%   combine instantaneous and window query points for a single kriging call;
%   hyperparameters are fit to the training data alone, so one call suffices ::
tqAll = [tq(:);        tqWindowFlat];
yqAll = [yq(:);        yqWindowFlat];
xqAll = [xq(:);        xqWindowFlat];

%   ordinary krige ::
[VhatEstMeanAnom, VhatStd] = ordinary_kriging(Xs(:, 2), ...
					      Xs(:, 1), ...
					      [], ...
					      anomV, ...
					      yqAll, ...
					      xqAll, ...
					      [], ...
					      'TrainTime', Xs(:, 3), ...
					      'QueryTime', tqAll, ...
					      'MaxK', 100, ...
					      'MinK', 20, ...
					      'FitHp', true, ...
					      'NRestarts', 3, ...
					      'RandomState', 7, ...
					      'Variance', 1.0, ...
					      'NoiseVariance', 0.1, ...
					      'GpHorizontalScaleM', 250000, ...
					      'GpTimeScale', 8.0, ...   % [d] MODIS 8-day composite period — physically motivated prior
					      'GpTimeScaleMin', 1.0, ...
					      'GpTimeScaleMax', 100.0);
% m2: old: 24.11/log(2) ≈ 34.8 d is Th-234 mean lifetime — copied from wrong context
% 'GpTimeScale', 24.11 / log(2), ...

%   invert temporal preprocessing for all combined query points ::
VhatEstMean = inverse_temporal_preprocess(VhatEstMeanAnom, tqAll, [], preproc);

%   split results: first nStations entries are instantaneous; remainder is window ::
VhatInstMean = VhatEstMean(1:nStations);
VhatInstStd  = VhatStd(1:nStations);

VhatWindowMean = VhatEstMean(nStations + 1 : end);
VhatWindowStd  = VhatStd(nStations + 1 : end);

%   reshape window results to nDaysWindow×nStations and average per station ::
%   column j of nppWindowMat contains the 35 daily predictions for station j
nppWindowMat    = reshape(VhatWindowMean, nDaysWindow, nStations);
nppStdWindowMat = reshape(VhatWindowStd,  nDaysWindow, nStations);

nppMean35d = mean(nppWindowMat, 1, 'omitnan')';   % nStations × 1
%   uncertainty: RMS of daily kriging SDs (conservative — positive temporal
%   correlation means the true mean SD is somewhat smaller, but this is safe)
nppKrigStd35d = sqrt(mean(nppStdWindowMat .^ 2, 1, 'omitnan'))';   % nStations × 1

%   store ::
nppStat = [sn, tq, VhatInstMean, VhatInstStd, nppMean35d, nppKrigStd35d];

%   make table ::
gp15_npp = array2table(nppStat);
gp15_npp.Properties.VariableNames = {'Station', 'Sampling_date', 'mmolC', 'stdev', 'mmolC_35d', 'stdev_35d'};

%   plot ::
figure; 
hold('on');
errorbar(gp15_stations.latitude, gp15_npp.mmolC, gp15_npp.stdev, '-ko', 'markerFaceColor', 'k', 'lineWidth', 1, 'markerSize', 10, 'capSize', 0); 
errorbar(gp15_stations.latitude, gp15_npp.mmolC_35d, gp15_npp.stdev_35d, '-ks', 'markerFaceColor', 'w', 'lineWidth', 1, 'markerSize', 10, 'capSize', 0); 
hold('off');
xlabel('\textbf{Latitude [deg N.]}', 'interpreter', 'latex');
ylabel('\textbf{NPP [mmol C m$^{-2}$ d$^{-1}$]}', 'interpreter', 'latex');
legend('Instantaneous', '35 [d] average', 'location', 'northWest', 'fontSize', 16, 'fontWeight', 'bold', 'interpreter', 'latex', 'box', 'off');
box('on');
set(gca, 'tickLabelInterpreter', 'latex', 'fontSize', 16, 'fontWeight', 'bold', 'lineWidth', 1);
title(['\textbf{CBPM NPP at GP15 Sampling Coordinates}'], 'interpreter', 'latex', 'fontSize', 20);
set(gcf, 'position', [0, 0, 1000, 500]); 
exportgraphics(gcf, [plot_output_basepath 'calcNpp/interp/gp15_npp.pdf'], 'ContentType', 'vector');

%%  make mean over cruise dates ::
%   get bounds ::
tMin = min(datenum(gp15_stations.date));
tMax = max(datenum(gp15_stations.date));
Ts = unique(T);

%   get sample bounds ::
[~, it1] = min(abs(tMax - Ts)); 
[~, it0] = min(abs(tMin - Ts)); 

%   make data ::
npp_mean = mean(NPP(:, :, it0:it1), 3, 'omitnan'); 
npp_se = std(NPP(:, :, it0:it1), 0, 3, 'omitnan') ./ sqrt(sum(isfinite(NPP(:, :, it0:it1)), 3));

%%  end subroutine
