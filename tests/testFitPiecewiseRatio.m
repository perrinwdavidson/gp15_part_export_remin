%%  testFitPiecewiseRatio - unit test for the piecewise ratio fitter
%   builds a synthetic dataset with known alpha, beta, z* and verifies
%   that fitPiecewiseRatio recovers them within tolerance.
%   the synthetic profile is exactly linear (noise added but within 5*sigma). ::

setupModel;

tol_ab   = 0.05;   % absolute tolerance on alpha, beta
tol_zstar = 50;    % absolute tolerance on z* [m]
nFail = 0;

%%  define ground truth
alpha_true = 5.0;
beta_true  = -0.02;
zstar_true = 200;
MLD        = 30;

%%  generate synthetic data: 20 stations x 5 depths
rng(42);
depths   = repmat([50; 100; 150; 200; 300], 20, 1);
stations = repelem((1:20)', 5);
noise    = 0.05 * randn(size(depths));
sigma    = 0.05 * ones(size(depths));

clampFn  = @(z) max(min(z, zstar_true), MLD);
ratio    = alpha_true + beta_true .* clampFn(depths) + noise;

%%  fit
% bug T7a: fBeta and fAlpha were swapped — function returns [alpha, beta, zstar, ...]
[fAlpha, fBeta, fZstar, ~, ~, flatLine, ~] = fitPiecewiseRatio(depths, ratio, sigma, stations, MLD);

%%  check
if abs(fBeta - beta_true) > tol_ab
    fprintf('  FAIL  beta: got %.4f, expected %.4f (tol %.4f)\n', fBeta, beta_true, tol_ab);
    nFail = nFail + 1;
else
    fprintf('  pass  beta (%.4f)\n', fBeta);
end

if abs(fAlpha - alpha_true) > tol_ab
    fprintf('  FAIL  alpha: got %.4f, expected %.4f (tol %.4f)\n', fAlpha, alpha_true, tol_ab);
    nFail = nFail + 1;
else
    fprintf('  pass  alpha (%.4f)\n', fAlpha);
end

if abs(fZstar - zstar_true) > tol_zstar
    fprintf('  FAIL  zstar: got %.1f m, expected %.1f m (tol %.1f m)\n', fZstar, zstar_true, tol_zstar);
    nFail = nFail + 1;
else
    fprintf('  pass  zstar (%.1f m)\n', fZstar);
end

if flatLine
    fprintf('  NOTE  flatLine=true (beta constrained to 0); check if noise is too large\n');
end

%%  report
if nFail == 0
    fprintf('testFitPiecewiseRatio: ALL PASS\n');
else
    fprintf('testFitPiecewiseRatio: %d FAILURE(S)\n', nFail);
end
