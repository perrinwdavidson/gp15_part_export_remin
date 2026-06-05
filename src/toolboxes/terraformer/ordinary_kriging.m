function [mean_pred, std_pred] = ordinary_kriging(train_lat, train_lon, train_depth, train_y, ...
                                                    query_lat, query_lon, query_depth, varargin)
% ORDINARY_KRIGING  True ordinary kriging via the Lagrange-multiplier system.
%
% Enforces the unbiasedness constraint sum(w) = 1 via a Lagrange multiplier.
% Returns the best linear unbiased predictor (BLUP) for an unknown constant
% mean field.  The returned std_pred is the true ordinary-kriging standard
% deviation, which is always >= the simple-kriging / GPR std.
%
% Augmented system solved for each query point x*:
%
%   [K(X,X)   1] [w ]   [k(x*, X)']
%   [1'       0] [lam] = [1        ]
%
% Distinction from simple_kriging
% --------------------------------
%   ordinary_kriging: weights sum to 1; var_OK >= var_SK; BLUP for unknown mean.
%   simple_kriging:   GP/SK with estimated mean; latent-field variance (slightly
%                     smaller std for small n).  Use for large-scale local GPR
%                     and future regression kriging (ML mean + SK residual).
%
% Syntax
% ------
%   [mean_pred, std_pred] = ordinary_kriging(train_lat, train_lon, train_depth, train_y, ...
%                                             query_lat, query_lon, query_depth, Name, Value)
%
% Coordinate semantics (mirrors Python _coerce_coord)
% ─────────────────────────────────────────────────────
%   Pass each spatial coordinate as:
%     numeric array  — active dimension included in the GP metric
%     scalar double  — constant reference (pinned value; not active in metric)
%     []             — excluded from metric entirely (Python None)
%
%   Supported dimension combinations include (but are not limited to):
%     Full spatial        : lat, lon, depth all arrays
%     Transect lat=const  : lat=scalar, lon+depth arrays
%     Transect lon=const  : lon=scalar, lat+depth arrays
%     Surface             : depth=[]
%     Lat+depth           : lon=[]
%     Depth-only          : lat=[], lon=[]
%     Mooring             : lat=[], lon=[], depth+time active
%     Time-only           : lat=[], lon=[], depth=[]  (requires TrainTime/QueryTime)
%
% Required positional arguments
% -----------------------------
%   train_lat, train_lon   coordinate (see semantics above).
%   train_depth            depth in metres, positive downward (array, scalar, or []).
%   train_y                [n×1 double]  Training observations (physical units).
%   query_lat, query_lon   query coordinates (same type as training coords).
%   query_depth            query depth (same type as train_depth).
%                          When a training coord is scalar, the query coord
%                          must be the same scalar value.
%
% Optional name–value pairs
% ─────────────────────────
% Spatiotemporal
%   'TrainTime'             [n×1 double]  Numeric training times in days.
%   'QueryTime'             [m×1 double]  Numeric query times in days.
%   'GpTimeScale'           scalar        GP temporal length scale (default 1095.75 days).
%   'GpTimeScaleMin'        scalar        Lower opt. bound (default 36.525 days).
%   'GpTimeScaleMax'        scalar        Upper opt. bound (default 36525 days).
%   'NeighborhoodTimeScale' scalar        BallTree temporal scale (default = GpTimeScale).
%
% Neighbourhood
%   'MaxK'         integer  Local mode: max neighbours per query.
%   'MinK'         integer  Min neighbours when radius is tight (default 3).
%   'SearchRadius' scalar   Optional radius filter (anisotropic ellipsoid semantics;
%                           see simple_kriging.m for full description).
%   'NHpFitMax'    integer  Max training pts for HP fitting (default 2000).
%
% Hyperparameters
%   'FitHp'         logical  Optimise hyperparameters (default true).
%   'Kernel'        string   'matern32'|'matern52'|'rbf' (default 'matern32').
%   'Variance'      scalar   GP signal variance (default 1.0).
%   'NoiseVariance' scalar   Observation noise variance (default 0.05).
%   'NRestarts'     integer  HP optimiser restarts (default 5).
%   'RandomState'   integer  Random seed (default = none).
%
% Spatial length scales
%   'GpHorizontalScaleM'  scalar (m)  Default 250000.
%   'GpVerticalScaleM'    scalar (m)  Default 100.
%
% Returns
% -------
%   mean_pred  [m×1]  BLUP for an unknown constant mean (weights sum to 1).
%   std_pred   [m×1]  True ordinary-kriging standard deviation (>= SK/GPR std).

p = inputParser();
p.CaseSensitive = false;

addParameter(p, 'TrainTime',   [],        @isnumeric);
addParameter(p, 'QueryTime',   [],        @isnumeric);
addParameter(p, 'GpTimeScale', 1095.75,   @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'GpTimeScaleMin', 36.525, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'GpTimeScaleMax', 36525., @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'NeighborhoodTimeScale', [], @isnumeric);
addParameter(p, 'MaxK',       [],     @(x) isempty(x) || (isnumeric(x) && isscalar(x)));
addParameter(p, 'MinK',       3,      @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'SearchRadius',[],    @(x) isempty(x) || (isnumeric(x) && isscalar(x)));
addParameter(p, 'NHpFitMax',  2000,   @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'FitHp',      true,   @islogical);
addParameter(p, 'Kernel',     'matern32', @ischar);
addParameter(p, 'Variance',   1.0,    @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'NoiseVariance', 0.05,@(x) isnumeric(x) && isscalar(x));
addParameter(p, 'NRestarts',  5,      @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'RandomState',[],     @(x) isempty(x) || isnumeric(x));
addParameter(p, 'GpHorizontalScaleM', 250000, @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'GpVerticalScaleM',   100,    @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'NeighborhoodHorizontalScaleM', [], @isnumeric);
addParameter(p, 'NeighborhoodVerticalScaleM',   [], @isnumeric);

parse(p, varargin{:});
opt = p.Results;

% ── Classify each coordinate axis (mirrors Python _coerce_coord) ──────────────
lat_is_scalar   = isnumeric(train_lat)   && isscalar(train_lat);
lat_is_none     = isempty(train_lat)     && ~isscalar(train_lat);
lat_is_array    = isnumeric(train_lat)   && ~isscalar(train_lat) && ~isempty(train_lat);

lon_is_scalar   = isnumeric(train_lon)   && isscalar(train_lon);
lon_is_none     = isempty(train_lon)     && ~isscalar(train_lon);
lon_is_array    = isnumeric(train_lon)   && ~isscalar(train_lon) && ~isempty(train_lon);

depth_is_scalar = isnumeric(train_depth) && isscalar(train_depth);
depth_is_none   = isempty(train_depth)   && ~isscalar(train_depth);
depth_is_array  = isnumeric(train_depth) && ~isscalar(train_depth) && ~isempty(train_depth);

% ── Validate query coordinate types match train semantics ─────────────────────
if lat_is_scalar
    if ~(isnumeric(query_lat) && isscalar(query_lat) && query_lat == train_lat)
        error('ordinary_kriging:queryLatMismatch', ...
              'query_lat must equal train_lat (%.15g) when lat is a scalar constant.', ...
              train_lat);
    end
elseif lat_is_none
    if ~isempty(query_lat)
        error('ordinary_kriging:queryLatMismatch', ...
              'query_lat must be [] when train_lat is [].');
    end
end

if lon_is_scalar
    if ~(isnumeric(query_lon) && isscalar(query_lon) && query_lon == train_lon)
        error('ordinary_kriging:queryLonMismatch', ...
              'query_lon must equal train_lon (%.15g) when lon is a scalar constant.', ...
              train_lon);
    end
elseif lon_is_none
    if ~isempty(query_lon)
        error('ordinary_kriging:queryLonMismatch', ...
              'query_lon must be [] when train_lon is [].');
    end
end

if depth_is_scalar
    if ~(isnumeric(query_depth) && isscalar(query_depth) && query_depth == train_depth)
        error('ordinary_kriging:queryDepthMismatch', ...
              'query_depth must equal train_depth (%.15g) when depth is a scalar constant.', ...
              train_depth);
    end
elseif depth_is_none
    if ~isempty(query_depth)
        error('ordinary_kriging:queryDepthMismatch', ...
              'query_depth must be [] when train_depth is [].');
    end
end

% ── Time ──────────────────────────────────────────────────────────────────────
has_train_time = ~isempty(opt.TrainTime);
has_query_time = ~isempty(opt.QueryTime);
if has_train_time ~= has_query_time
    error('ordinary_kriging:timeMismatch', ...
          'Both TrainTime and QueryTime must be provided, or neither.');
end

if isempty(opt.NeighborhoodTimeScale)
    nt_scale = opt.GpTimeScale;
else
    nt_scale = opt.NeighborhoodTimeScale;
end

% ── Write CSVs ─────────────────────────────────────────────────────────────────
train_file  = [tempname, '.csv'];
query_file  = [tempname, '.csv'];
result_file = [tempname, '.csv'];
cleanup = onCleanup(@() delete_if_exists({train_file, query_file, result_file}));

train_var_names = {'y'};
train_arrays    = {train_y(:)};
if lat_is_array
    train_var_names = [train_var_names, {'lat'}];
    train_arrays    = [train_arrays,    {train_lat(:)}];
end
if lon_is_array
    train_var_names = [train_var_names, {'lon'}];
    train_arrays    = [train_arrays,    {train_lon(:)}];
end
if depth_is_array
    train_var_names = [train_var_names, {'depth'}];
    train_arrays    = [train_arrays,    {train_depth(:)}];
end
if has_train_time
    train_var_names = [train_var_names, {'time'}];
    train_arrays    = [train_arrays,    {opt.TrainTime(:)}];
end
T_train = table(train_arrays{:}, 'VariableNames', train_var_names);

query_var_names = {};
query_arrays    = {};
if lat_is_array
    query_var_names = [query_var_names, {'lat'}];
    query_arrays    = [query_arrays,    {query_lat(:)}];
end
if lon_is_array
    query_var_names = [query_var_names, {'lon'}];
    query_arrays    = [query_arrays,    {query_lon(:)}];
end
if depth_is_array
    query_var_names = [query_var_names, {'depth'}];
    query_arrays    = [query_arrays,    {query_depth(:)}];
end
if has_train_time
    query_var_names = [query_var_names, {'time'}];
    query_arrays    = [query_arrays,    {opt.QueryTime(:)}];
end

if isempty(query_var_names)
    error('ordinary_kriging:noActiveDims', ...
          ['No active dimensions in query: all coordinates are scalar/[] ' ...
           'and no time was provided.  At least one active dimension is required.']);
end
T_query = table(query_arrays{:}, 'VariableNames', query_var_names);

writetable(T_train, train_file);
writetable(T_query, query_file);

% ── Build CLI command (calls terraformer._cli.ok = true OK engine) ────────────
cmd = sprintf('python -m terraformer._cli.ok --train_csv "%s" --query_csv "%s" --result_csv "%s"', ...
              train_file, query_file, result_file);

if lat_is_scalar
    cmd = [cmd, sprintf(' --lat_scalar %.15g', train_lat)];
elseif lat_is_none
    cmd = [cmd, ' --no_lat'];
end
if lon_is_scalar
    cmd = [cmd, sprintf(' --lon_scalar %.15g', train_lon)];
elseif lon_is_none
    cmd = [cmd, ' --no_lon'];
end
if depth_is_scalar
    cmd = [cmd, sprintf(' --depth_scalar %.15g', train_depth)];
elseif depth_is_none
    cmd = [cmd, ' --no_depth'];
end

if ~isempty(opt.MaxK)
    cmd = [cmd, sprintf(' --max_k %d', int32(opt.MaxK))];
end
cmd = [cmd, sprintf(' --min_k %d --n_hp_fit_max %d', int32(opt.MinK), int32(opt.NHpFitMax))];
if ~isempty(opt.SearchRadius)
    cmd = [cmd, sprintf(' --search_radius %.6g', opt.SearchRadius)];
end

if ~opt.FitHp,  cmd = [cmd, ' --no_fit_hp'];  end
cmd = [cmd, sprintf(' --kernel %s --variance %.10g --noise_variance %.10g', ...
                    opt.Kernel, opt.Variance, opt.NoiseVariance)];
cmd = [cmd, sprintf(' --n_restarts %d', int32(opt.NRestarts))];
if ~isempty(opt.RandomState)
    cmd = [cmd, sprintf(' --random_state %d', int32(opt.RandomState))];
end

cmd = [cmd, sprintf(' --gp_horizontal_scale_m %.6g --gp_vertical_scale_m %.6g', ...
                    opt.GpHorizontalScaleM, opt.GpVerticalScaleM)];
if ~isempty(opt.NeighborhoodHorizontalScaleM)
    cmd = [cmd, sprintf(' --neighborhood_horizontal_scale_m %.6g', opt.NeighborhoodHorizontalScaleM)];
end
if ~isempty(opt.NeighborhoodVerticalScaleM)
    cmd = [cmd, sprintf(' --neighborhood_vertical_scale_m %.6g', opt.NeighborhoodVerticalScaleM)];
end
cmd = [cmd, sprintf(' --gp_time_scale %.10g --gp_time_scale_min %.10g --gp_time_scale_max %.10g', ...
                    opt.GpTimeScale, opt.GpTimeScaleMin, opt.GpTimeScaleMax)];
cmd = [cmd, sprintf(' --neighborhood_time_scale %.10g', nt_scale)];

[status, cmdout] = system(cmd);
if status ~= 0
    error('ordinary_kriging:pythonError', ...
          'terraformer._cli.ok returned status %d:\n%s', status, cmdout);
end
if ~isfile(result_file)
    error('ordinary_kriging:noResultFile', ...
          'Result CSV not created.\n%s', cmdout);
end

T_result  = readtable(result_file);
mean_pred = T_result.mean(:);
std_pred  = T_result.std(:);
end

function delete_if_exists(files)
    for i = 1:numel(files)
        if isfile(files{i}),  delete(files{i});  end
    end
end
