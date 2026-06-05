function mean_physical = inverse_temporal_preprocess(mean_anomaly, query_time, query_depth, preproc_path)
% INVERSE_TEMPORAL_PREPROCESS  Convert GP anomaly predictions to physical units.
%
% Syntax
% ------
%   mean_physical = inverse_temporal_preprocess(mean_anomaly, query_time, ...
%                                               query_depth, preproc_path)
%
% Description
% -----------
%   Applies the inverse of a fitted TemporalPreprocessor to GP-predicted
%   anomaly values, adding back the seasonal+trend component to recover
%   physical units.
%
%   This is the companion to fit_temporal_preprocessor — pass the
%   preproc_path returned by that function as the final argument.
%
% Arguments
% ---------
%   mean_anomaly  [m×1 double]  GP posterior mean in anomaly space.
%   query_time    [m×1 double]  Query times (float64 days, same epoch as training).
%   query_depth   [m×1 double]  Query depths (metres, positive downward).
%                               Pass [] to use the depth-free global model.
%   preproc_path  char          Path to .joblib file from fit_temporal_preprocessor.
%
% Returns
% -------
%   mean_physical  [m×1 double]  GP posterior mean in physical units.

if nargin < 3
    query_depth = [];
end

% ── Paths ─────────────────────────────────────────────────────────────────────
gp_csv     = [tempname, '.csv'];
query_csv  = [tempname, '.csv'];
output_csv = [tempname, '.csv'];
cleanup = onCleanup(@() itp_delete_temps({gp_csv, query_csv, output_csv}));

% ── Write GP output CSV (mean + zero std placeholder) ─────────────────────────
T_gp = table(mean_anomaly(:), zeros(numel(mean_anomaly), 1), ...
             'VariableNames', {'mean', 'std'});
writetable(T_gp, gp_csv);

% ── Write query CSV ───────────────────────────────────────────────────────────
if isempty(query_depth)
    T_q = table(query_time(:), 'VariableNames', {'time'});
else
    T_q = table(query_time(:), query_depth(:), 'VariableNames', {'time', 'depth'});
end
writetable(T_q, query_csv);

% ── Build CLI command ─────────────────────────────────────────────────────────
cmd = sprintf( ...
    'python -m terraformer._cli.preproc --mode inverse --preprocessor_path "%s" --input_csv "%s" --query_csv "%s" --output_csv "%s"', ...
    preproc_path, gp_csv, query_csv, output_csv);

% ── Run ───────────────────────────────────────────────────────────────────────
[status, cmdout] = system(cmd);
if status ~= 0
    error('inverse_temporal_preprocess:pythonError', ...
          'terraformer._cli.preproc returned status %d:\n%s', status, cmdout);
end
if ~isfile(output_csv)
    error('inverse_temporal_preprocess:noOutputFile', ...
          'Output CSV not created.\n%s', cmdout);
end

% ── Read result ───────────────────────────────────────────────────────────────
T_out         = readtable(output_csv);
mean_physical = T_out.mean(:);
end

function itp_delete_temps(files)
    for i = 1:numel(files)
        if isfile(files{i}),  delete(files{i});  end
    end
end
