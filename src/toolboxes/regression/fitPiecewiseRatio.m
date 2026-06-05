function [fittedAlpha, fittedBeta, fittedZStar, paramCov, zStarCI, flatLine, nBinMeans] = ...
         fitPiecewiseRatio(z_raw, r_raw, sigma_raw, station_raw, mld)
%% fitPiecewiseRatio  profile-likelihood piecewise ratio regression
%
%  Fits the three-zone model with a hard non-positive slope constraint:
%
%    R(z) = alpha + beta * clamp(z, mld, zStar),   beta <= 0
%
%  where clamp(z, lo, hi) = max(min(z, hi), lo).  The model is constant
%  above mld (surface cap), linearly decreasing from mld to zStar, then
%  constant below zStar — continuous at every boundary by construction.
%
%  The transition depth zStar is found by sweeping a profile likelihood
%  over the unique raw measurement depths in [mld, BTMDEPTH].  At each
%  candidate zStar the inner (alpha, beta) problem is a constrained WLS
%  that is solvable in closed form: accept the unconstrained solution when
%  beta_unc <= 0, otherwise clamp to the boundary beta = 0 (flat line,
%  alpha = weighted mean).
%
%  To reduce within-station serial correlation, per-station inverse-
%  variance weighted means are computed within equal-count depth bins
%  before fitting.  The number of bins is chosen as
%    n_bins = max(5, min(10, floor(n_raw / 5)))
%  so each bin holds at least ~5 raw observations.
%
%  Uncertainty on (alpha, beta) uses the WLS covariance scaled by the
%  residual variance (sSquared = RSS / dof), consistent with ols.m.
%  The 95% CI on zStar comes from the concentrated profile likelihood:
%    CI = { zs : n * log(RSS(zs)/RSS_opt) <= chi2(0.95,1) = 3.841 }
%
%  Inputs
%    z_raw       - depths [m]                   (n_raw x 1)
%    r_raw       - ratio values                 (n_raw x 1)
%    sigma_raw   - 1-sigma ratio uncertainties  (n_raw x 1)
%    station_raw - station label per obs        (n_raw x 1)
%    mld         - regional mean MLD [m]; lower bound for zStar search
%
%  Outputs
%    fittedAlpha  - intercept  [same units as r_raw]
%    fittedBeta   - slope [units/m],  <= 0  (== 0 when flatLine is true)
%    fittedZStar  - transition depth [m]
%    paramCov     - 2x2 covariance of [alpha; beta] at fittedZStar
%    zStarCI      - [lo, hi] 95% profile-likelihood CI on zStar [m]
%    flatLine     - true when the beta <= 0 constraint is active at zero
%    nBinMeans    - number of per-station-per-bin means used in the fit

    BTMDEPTH = 400;    % [m] — must equal the pipeline BTMDEPTH constant
    CHI2_95  = 3.841;  % chi-squared 95th percentile, 1 dof

    z_raw       = z_raw(:);
    r_raw       = r_raw(:);
    sigma_raw   = sigma_raw(:);
    station_raw = station_raw(:);

    % ----------------------------------------------------------------
    % 1.  Equal-count depth binning
    % ----------------------------------------------------------------
    n_raw  = length(z_raw);
    n_bins = max(5, min(10, floor(n_raw / 5)));

    % Rank-based assignment: each bin receives floor(n_raw/n_bins)
    % observations regardless of depth ties
    [~, sort_ord]     = sort(z_raw);
    bin_of_rank       = min(ceil((1:n_raw)' * n_bins / n_raw), n_bins);
    bin_idx           = zeros(n_raw, 1);
    bin_idx(sort_ord) = bin_of_rank;

    % ----------------------------------------------------------------
    % 2.  Per-station-per-bin inverse-variance weighted means
    % ----------------------------------------------------------------
    uniqueStats = unique(station_raw);
    z_bin   = zeros(0, 1);
    r_bin   = zeros(0, 1);
    sig_bin = zeros(0, 1);

    for s = 1:length(uniqueStats)
        sIdx = (station_raw == uniqueStats(s));
        for b = 1:n_bins
            sbIdx = sIdx & (bin_idx == b);
            if ~any(sbIdx), continue; end

            z_sb  = z_raw(sbIdx);
            r_sb  = r_raw(sbIdx);
            sg_sb = sigma_raw(sbIdx);

            ok = isfinite(r_sb) & isfinite(sg_sb) & (sg_sb > 0);
            if ~any(ok), continue; end
            z_sb  = z_sb(ok);
            r_sb  = r_sb(ok);
            sg_sb = sg_sb(ok);

            w_sb    = 1 ./ sg_sb.^2;
            z_bin   = [z_bin;   mean(z_sb)];
            r_bin   = [r_bin;   sum(w_sb .* r_sb) / sum(w_sb)];
            sig_bin = [sig_bin; 1 / sqrt(sum(w_sb))];
        end
    end

    nBinMeans = length(z_bin);
    if nBinMeans < 3
        error('fitPiecewiseRatio: only %d bin means — need at least 3.', nBinMeans);
    end
    w_bin = 1 ./ sig_bin.^2;

    % ----------------------------------------------------------------
    % 3.  Profile likelihood: sweep over candidate zStar values
    % ----------------------------------------------------------------
    z_cands = unique(z_raw(z_raw >= mld & z_raw <= BTMDEPTH));
    if isempty(z_cands)
        error('fitPiecewiseRatio: no raw measurements in [mld=%.1f, %d m].', mld, BTMDEPTH);
    end
    nCands  = length(z_cands);
    RSS_vec = zeros(nCands, 1);
    a_vec   = zeros(nCands, 1);
    b_vec   = zeros(nCands, 1);

    for k = 1:nCands
        [a_vec(k), b_vec(k), RSS_vec(k)] = innerWLS(z_bin, r_bin, w_bin, z_cands(k), mld);
    end

    % ----------------------------------------------------------------
    % 4.  Optimal parameters and covariance
    % ----------------------------------------------------------------
    [RSS_opt, k_opt] = min(RSS_vec);
    fittedZStar = z_cands(k_opt);
    fittedAlpha = a_vec(k_opt);
    fittedBeta  = b_vec(k_opt);
    flatLine    = (fittedBeta == 0);

    zActive_opt = max(min(z_bin, fittedZStar), mld);
    X_opt       = [ones(nBinMeans, 1), zActive_opt];
    if ~flatLine
        sSquared = RSS_opt / max(nBinMeans - 2, 1);
        XtWX     = X_opt' * (X_opt .* w_bin);
        paramCov = sSquared * inv(XtWX);
    else
        % beta fixed at zero boundary: effective 1-parameter model
        sSquared = RSS_opt / max(nBinMeans - 1, 1);
        paramCov = [sSquared / sum(w_bin), 0; 0, 0];
    end

    % ----------------------------------------------------------------
    % 5.  95% profile-likelihood CI on zStar
    %     Concentrated likelihood: -2 logLR = n * log(RSS(zs)/RSS_opt)
    %     CI set: RSS(zs) <= RSS_opt * exp(3.841 / n)
    % ----------------------------------------------------------------
    ci_thresh = RSS_opt * exp(CHI2_95 / nBinMeans);
    in_ci     = (RSS_vec <= ci_thresh);
    if any(in_ci)
        zStarCI = [min(z_cands(in_ci)), max(z_cands(in_ci))];
    else
        zStarCI = [fittedZStar, fittedZStar];
    end

end

% -----------------------------------------------------------------------
function [alpha, beta, RSS] = innerWLS(z_bin, r_bin, w_bin, zStar, mld)
%  Constrained WLS for fixed zStar.
%  Minimises sum_i w_i * (r_i - alpha - beta*clamp(z_i,mld,zStar))^2
%  subject to beta <= 0.  Solved analytically: the feasible set {beta<=0}
%  is a half-space, so the constrained optimum is either the unconstrained
%  solution (when beta_unc <= 0) or the boundary solution beta=0.

    zActive = max(min(z_bin, zStar), mld);
    X       = [ones(length(z_bin), 1), zActive];

    Xw    = X .* sqrt(w_bin);
    rw    = r_bin .* sqrt(w_bin);
    theta = Xw \ rw;

    if theta(2) <= 0
        alpha = theta(1);
        beta  = theta(2);
    else
        % constraint active: flat line, alpha = weighted mean
        beta  = 0;
        alpha = sum(w_bin .* r_bin) / sum(w_bin);
    end

    RSS = sum(w_bin .* (r_bin - alpha - beta * zActive).^2);
end
