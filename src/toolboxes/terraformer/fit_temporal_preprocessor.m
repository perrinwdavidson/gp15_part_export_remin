function [anomaly_y, preproc_path] = fit_temporal_preprocessor(lat, lon, depth, y, time, varargin)
% FIT_TEMPORAL_PREPROCESSOR  Fit a TemporalPreprocessor and return anomaly y.
%
% Syntax
% ------
%   [anomaly_y, preproc_path] = fit_temporal_preprocessor([], [], [], y, time)
%   [anomaly_y, preproc_path] = fit_temporal_preprocessor(lat, lon, depth, y, time)
%   [anomaly_y, preproc_path] = fit_temporal_preprocessor(lat, lon, depth, y, time, ...
%                                   'DepthBinEdges', edges, 'MinSamplesPerBin', 30)
%   [anomaly_y, preproc_path] = fit_temporal_preprocessor(lat, lon, depth, y, time, ...
%                                   'MaxSamplesPerBin', 500, 'RandomState', 42)
%
% Description
% -----------
%   Fits a depth-binned harmonic regression (annual + semi-annual cycles +
%   linear trend) on the training observations, subtracts the fitted trend to
%   produce anomalies, and saves the fitted preprocessor to a .joblib file for
%   later use with inverse_temporal_preprocess.
%
%   Preprocessing is intentionally separate from simple_kriging / ordinary_kriging
%   — pass anomaly_y to the kriging function and use inverse_temporal_preprocess
%   to convert GP predictions back to physical units.
%
%   Pipeline:
%       [anom_y, pp] = fit_temporal_preprocessor(lat, lon, dep, y, t);
%       [mean_anom, std] = simple_kriging(lat, lon, dep, anom_y, ...
%                              q_lat, q_lon, q_dep, 'TrainTime', t, ...
%                              'QueryTime', q_t, ...);
%       mean_phys = inverse_temporal_preprocess(mean_anom, q_t, q_dep, pp);
%
% Required arguments
% ------------------
%   y      [n×1 double]   Target observations (physical units).
%   time   [n×1 double]   Numeric times in float64 days since epoch.
%
% Optional positional arguments (pass [] to omit)
% ------------------------------------------------
%   lat    [n×1 double] or []   Latitudes (degrees).  Not used by the
%                               preprocessor; pass [] when unavailable.
%   lon    [n×1 double] or []   Longitudes (degrees).  Same as lat.
%   depth  [n×1 double] or []   Depths (metres, positive downward).
%                               When [] or omitted: one global harmonic model
%                               is fitted across all observations (no depth
%                               binning).  When a constant scalar is given
%                               (e.g., all data at 200 m), the bin containing
%                               that depth is used.
%
% Optional name–value pairs
% ─────────────────────────
%   'DepthBinEdges'     [1×k double]  Depth bin edges (metres).
%                                     Default: [0, 200, 1000, Inf].
%   'MinSamplesPerBin'  integer       Min samples per bin before using the
%                                     global model as fallback. Default: 30.
%   'MaxSamplesPerBin'  integer or [] Cap on the number of observations used
%                                     for each OLS fit (both the global model
%                                     and each per-bin model).  When a bin or
%                                     the global dataset exceeds this count, a
%                                     random subsample of this size is drawn
%                                     before fitting.  Use with 'RandomState'
%                                     for reproducibility.  Default: [] (use
%                                     all observations).
%   'RandomState'       integer or [] Integer seed for the subsampling RNG.
%                                     Only has an effect when 'MaxSamplesPerBin'
%                                     is set.  Default: [] (non-deterministic).
%   'PreprocessorPath'  char          Destination for the saved .joblib.
%                                     Default: auto-generated tempname.
%
% Returns
% -------
%   anomaly_y    [n×1 double]  Seasonal-trend-removed observations.
%   preproc_path char          Path to the saved .joblib file.

p = inputParser();
p.CaseSensitive = false;
addParameter(p, 'DepthBinEdges',    [],  @isnumeric);
addParameter(p, 'MinSamplesPerBin', [],  @(x) isempty(x) || (isnumeric(x) && isscalar(x)));
addParameter(p, 'MaxSamplesPerBin', [],  @(x) isempty(x) || (isnumeric(x) && isscalar(x)));
addParameter(p, 'RandomState',      [],  @(x) isempty(x) || (isnumeric(x) && isscalar(x)));
addParameter(p, 'PreprocessorPath', '',  @ischar);
parse(p, varargin{:});
opt = p.Results;

% ── Paths ─────────────────────────────────────────────────────────────────────
if isempty(opt.PreprocessorPath)
    preproc_path = [tempname, '.joblib'];
else
    preproc_path = opt.PreprocessorPath;
end
input_csv   = [tempname, '.csv'];
anomaly_csv = [tempname, '.csv'];
cleanup = onCleanup(@() ftp_delete_temps({input_csv, anomaly_csv}));

% ── Write input CSV ───────────────────────────────────────────────────────────
% y and time are always written; lat, lon, depth are included only when
% non-empty so the Python CLI can detect which dimensions are available.
T = table(y(:), time(:), 'VariableNames', {'y', 'time'});
if ~isempty(lat),   T.lat   = lat(:);   end
if ~isempty(lon),   T.lon   = lon(:);   end
if ~isempty(depth), T.depth = depth(:); end
writetable(T, input_csv);

% ── Build CLI command ─────────────────────────────────────────────────────────
cmd = sprintf('python -m terraformer._cli.preproc --mode fit --input_csv "%s" --anomaly_csv "%s" --preprocessor_path "%s"', ...
              input_csv, anomaly_csv, preproc_path);

if ~isempty(opt.DepthBinEdges)
    edges_str = strjoin(arrayfun(@(x) num2str(x, '%.6g'), opt.DepthBinEdges(:)', ...
                                 'UniformOutput', false), ',');
    edges_str = strrep(edges_str, 'Inf', 'inf');
    cmd = [cmd, sprintf(' --depth_bin_edges "%s"', edges_str)];
end
if ~isempty(opt.MinSamplesPerBin)
    cmd = [cmd, sprintf(' --min_samples_per_bin %d', int32(opt.MinSamplesPerBin))];
end
if ~isempty(opt.MaxSamplesPerBin)
    cmd = [cmd, sprintf(' --max_samples_per_bin %d', int32(opt.MaxSamplesPerBin))];
end
if ~isempty(opt.RandomState)
    cmd = [cmd, sprintf(' --random_state %d', int32(opt.RandomState))];
end

% ── Run ───────────────────────────────────────────────────────────────────────
[status, cmdout] = system(cmd);
if status ~= 0
    error('fit_temporal_preprocessor:pythonError', ...
          'terraformer._cli.preproc returned status %d:\n%s', status, cmdout);
end
if ~isfile(anomaly_csv)
    error('fit_temporal_preprocessor:noOutputFile', ...
          'Anomaly CSV not created.\n%s', cmdout);
end

% ── Read anomaly CSV ──────────────────────────────────────────────────────────
T_anom    = readtable(anomaly_csv);
anomaly_y = T_anom.y(:);
end

function ftp_delete_temps(files)
    for i = 1:numel(files)
        if isfile(files{i}),  delete(files{i});  end
    end
end
