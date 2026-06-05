%%  testCalcError - unit test for the 1D and 2D flux variance formulae
%   hand-computes expected variances for a 3-depth synthetic station and
%   compares against the pipeline formula using calcLayerThickness.
%   pass criteria: relative error < 1e-12. ::

setupModel;

tol = 1e-12;
nFail = 0;

%%  define synthetic 3-depth station
z   = [10; 50; 100];       % depths [m]
nu  = [0.1; 0.1; 0.1];    % sigma_U238 [dpm m-3] (in model units)
nth = [0.2; 0.2; 0.2];    % sigma_Th234 [dpm m-3]
w   = [0.001; 0.002; 0.003];     % upwelling velocity [m d-1]
nw  = [1e-4;  1e-4;  1e-4];     % sigma_w [m d-1]
g   = [0.010; 0.020; 0.015];    % vertical gradient [dpm m-4]
ng  = [0.001; 0.001; 0.001];    % sigma_gradient [dpm m-4]

%%  compute dz via calcLayerThickness
dz = calcLayerThickness(z);
% expected: dz = [30; 45; 25] for z = [10; 50; 100]
dz_expected = [30; 45; 25];
nFail = nFail + checkVal('dz', dz, dz_expected, tol);

%%  compute expected 1D variance (hand calculation)
LAMBDA = log(2) / 24.101;
perLayer1d = (LAMBDA .* dz) .^ 2 .* (nu .^ 2 + nth .^ 2);
error1d_expected = cumsum(perLayer1d);

%   formula from calcError.m ::
error1d = cumsum((((LAMBDA .* dz) .^ 2) .* ((nu .^ 2) + (nth .^ 2))));
nFail = nFail + checkVal('error1d', error1d, error1d_expected, tol);

%%  compute expected 2D variance (hand calculation)
upwellError_expected = (g .* dz) .^ 2 .* nw .^ 2 + (w .* dz) .^ 2 .* ng .^ 2;
error2d_expected     = cumsum(perLayer1d + upwellError_expected);

%   formula from calcError.m (absolute-error form, E1 fix) ::
upwellError = (g .* dz) .^ 2 .* nw .^ 2 + (w .* dz) .^ 2 .* ng .^ 2;
error2d     = cumsum((((LAMBDA .* dz) .^ 2) .* ((nu .^ 2) + (nth .^ 2))) + upwellError);
nFail = nFail + checkVal('error2d', error2d, error2d_expected, tol);

%%  report
if nFail == 0
    fprintf('testCalcError: ALL PASS\n');
else
    fprintf('testCalcError: %d FAILURE(S)\n', nFail);
end

%%  helper
function fail = checkVal(name, actual, expected, tol)
    relerr = max(abs(actual - expected) ./ (abs(expected) + eps));
    if relerr > tol
        fprintf('  FAIL  %s: max relative error = %.2e (tol = %.2e)\n', name, relerr, tol);
        fail = 1;
    else
        fprintf('  pass  %s\n', name);
        fail = 0;
    end
end
