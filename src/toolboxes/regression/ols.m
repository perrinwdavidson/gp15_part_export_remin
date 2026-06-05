%% ols -- ordinary / weighted least squares regression
%% --------------------------------------------------------
function [beta, alpha, q, ci] = ols(x, y, y_uncertainty, plot_data, test_data)

    % Defaults
    if nargin < 3
        y_uncertainty = [];
    end

    if nargin < 4 || isempty(plot_data)
        plot_data = false;
    end

    if nargin < 5 || isempty(test_data)
        test_data = false;
    end

    % Optional test data
    if test_data
        sz = [1E2, 1];
        x = rand(sz);
        epsilon = normrnd(0, 1E-1, sz);
        y = x + epsilon;

        if isempty(y_uncertainty)
            y_uncertainty = 1E-1 * ones(size(y));
        end
    end

    % Force column vectors
    x = x(:);
    y = y(:);

    % Check lengths
    if length(x) ~= length(y)
        error('Input vectors x and y are not the same length.');
    end

    % Determine whether to use WLS
    use_wls = ~isempty(y_uncertainty);

    if use_wls
        y_uncertainty = y_uncertainty(:);

        if length(y_uncertainty) ~= length(y)
            error('y_uncertainty must have the same length as y.');
        end

        if any(y_uncertainty <= 0 | ~isfinite(y_uncertainty))
            error('All y_uncertainty values must be positive and finite.');
        end
    end

    % Remove NaNs / Infs
    if use_wls
        valid = isfinite(x) & isfinite(y) & isfinite(y_uncertainty);
        x = x(valid);
        y = y(valid);
        y_uncertainty = y_uncertainty(valid);
    else
        valid = isfinite(x) & isfinite(y);
        x = x(valid);
        y = y(valid);
    end

    % Number of data points
    n = length(x);
    dof = n - 2;

    if n < 3
        error('Need at least 3 finite data points for regression.');
    end

    % Design matrix
    % y = alpha + beta*x
    X = [ones(n, 1), x];

    % OLS or WLS fit
    if use_wls
        w = 1 ./ (y_uncertainty .^ 2);

        % Avoid explicitly forming W = diag(w)
        Xw = X .* sqrt(w);
        yw = y .* sqrt(w);

        theta = Xw \ yw;

        y_fit = X * theta;
        RSS = sum(w .* (y - y_fit).^2);

        sSquared = RSS / dof;

        covTheta = sSquared * inv(X' * (X .* w));

        fit_type = 'WLS';
    else
        theta = X \ y;

        y_fit = X * theta;
        RSS = sum((y - y_fit).^2);

        sSquared = RSS / dof;

        covTheta = sSquared * inv(X' * X);

        fit_type = 'OLS';
    end

    % Extract intercept and slope
    alpha = theta(1);
    beta  = theta(2);

    % Standard error of fitted mean at each observed x
    q = zeros(size(x));

    for i = 1:n
        xi = [1, x(i)];
        q(i) = sqrt(xi * covTheta * xi');
    end

    % Store useful regression diagnostics
    Ci = struct();
    Ci.varAlpha = covTheta(1, 1);
    Ci.varBeta  = covTheta(2, 2);
    Ci.covAlphaBeta = covTheta(1, 2);
    Ci.n = n;
    Ci.dof = dof;
    Ci.RSS = RSS;
    Ci.sSquared = sSquared;
    Ci.covTheta = covTheta;
    Ci.fit_type = fit_type;
    Ci.use_wls = use_wls;
    % ci = [varBeta, varAlpha, covAlphaBeta, n];  % bug: bare names undefined; values live in Ci struct
    ci = [Ci.varBeta, Ci.varAlpha, Ci.covAlphaBeta, n];

    % Plot
    if plot_data
        figure;
        hold on;

        if use_wls
            errorbar(x, y, y_uncertainty, 'ko', ...
                'LineWidth', 1.5, ...
                'MarkerFaceColor', 'k');
        else
            scatter(x, y, 80, 'k', 'filled');
        end

        xlabel('x');
        ylabel('y');
        title(sprintf('%s Regression', fit_type));

        x_tilde = linspace(min(x), max(x), 1E3)';
        X_tilde = [ones(length(x_tilde), 1), x_tilde];

        y_tilde = X_tilde * theta;

        q_tilde = zeros(size(x_tilde));

        for i = 1:length(x_tilde)
            xi = X_tilde(i, :);
            q_tilde(i) = sqrt(xi * covTheta * xi');
        end

        y_lower = y_tilde - q_tilde;
        y_upper = y_tilde + q_tilde;

        plot(x_tilde, y_tilde, 'r', 'LineWidth', 2);

        fill([x_tilde; flipud(x_tilde)], ...
             [y_lower; flipud(y_upper)], ...
             'r', ...
             'FaceAlpha', 0.2, ...
             'EdgeColor', 'none');

        hold off;

        if use_wls
            legend('Data with uncertainty', fit_type, 'SE', ...
                   'Location', 'southeast');
        else
            legend('Data', fit_type, 'SE', ...
                   'Location', 'southeast');
        end

        set(gca, 'Box', 'on', ...
                 'XLim', [min(x), max(x)]);
    end

end
