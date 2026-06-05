%%  testUncertMetrics - unit test for T100 and R100 uncertainty formulae
%   Synthetic 2-layer station with known per-layer variances. Verifies:
%     R100: sigma = sqrt(Var[incremental layer])  (nested covariance)
%     T100: sigma from Th-234 variances — exact because ratio R cancels algebraically
%   All expected values are hand-derived from first principles; pass criteria: rel-err < 1e-12. ::

setupModel;

tol   = 1e-12;
nFail = 0;

%%  synthetic 2-layer station
%   Layer 1 (surface to PPZ):   Th contribution f1 = 100, Var[f1] = 100
%   Layer 2 (PPZ to PPZ+100):   Th contribution f2 = -30, Var[f2] = 36
%                                (negative = net Th excess from remineralisation below PPZ)
%   Deep-zone ratio (same constant at both depths): R = 0.04, sigma_R/R = 0.15

% cumulative Th-234 fluxes
F_Th_PPZ    = 100;          % dpm m-2 d-1
F_Th_100Ppz = 70;           % = f1 + f2

% cumulative Th-234 variances (independent layers -> additive)
Var_Th_PPZ    = 100;        % = Var[f1]
Var_Th_100Ppz = 100 + 36;   % = Var[f1] + Var[f2]   ->  136

% ratio and its uncertainty (same constant at both depths in deep zone)
R         = 0.04;           % mmol dpm-1
relSigmaR = 0.15;           % sigma_R / R
sigmaR    = R * relSigmaR;  % 0.006

% POC fluxes
F_POC_PPZ    = F_Th_PPZ    * R;   % 4.0 mmol m-2 d-1
F_POC_100Ppz = F_Th_100Ppz * R;   % 2.8 mmol m-2 d-1

% POC variances via delta method for product F_POC = F_Th * R (F_Th, R independent)
Var_POC_PPZ    = R^2 * Var_Th_PPZ    + F_Th_PPZ^2    * sigmaR^2;   % 0.16 + 0.36  = 0.52
Var_POC_100Ppz = R^2 * Var_Th_100Ppz + F_Th_100Ppz^2 * sigmaR^2;   % 0.2176+0.1764 = 0.394

%%  R100: nested covariance for Th-234 cumulative sums
%   Cov[F_Th(PPZ+100), F_Th(PPZ)] = Var[F_Th(PPZ)]  (shared surface-to-PPZ layers)
%   => Var[R100] = Var[F_Th(PPZ+100)] - Var[F_Th(PPZ)] = Var[f2] = 36
%   => sigma_R100 = 6

R100 = F_Th_100Ppz - F_Th_PPZ;
nFail = nFail + checkVal('R100 central value',  R100,                                   -30, tol);
nFail = nFail + checkVal('R100 uncertainty',    sqrt(max(Var_Th_100Ppz - Var_Th_PPZ, 0)), 6, tol);

%%  T100: POC central value, Th-234 variance formula
%   T100 = F_POC(PPZ+100)/F_POC(PPZ) = 2.8/4.0 = 0.7
%   In the delta method for T100 = F_Th(PPZ+100)*R / (F_Th(PPZ)*R), with R constant at
%   both depths and independent of F_Th, the exact covariance is:
%     Cov[a,b] = R^2*Var_Th(PPZ) + F_Th(PPZ+100)*F_Th(PPZ)*sigma_R^2
%   Substituting into (Var[a] + T100^2*Var[b] - 2*T100*Cov)/b^2, the sigma_R^2 terms
%   cancel exactly, giving:
%     Var[T100] = (Var_Th(PPZ+100) - T100*(2-T100)*Var_Th(PPZ)) / F_Th(PPZ)^2
%              = (136 - 0.7*1.3*100) / 100^2 = 45/10000 = 0.0045
%   sigma_T100 = sqrt(0.0045) = sqrt(45)/100

T100 = F_POC_100Ppz / F_POC_PPZ;
nFail = nFail + checkVal('T100 central value', T100, 0.7, tol);

sigma_Th_PPZ    = sqrt(Var_Th_PPZ);
sigma_Th_100Ppz = sqrt(Var_Th_100Ppz);
uncertT100      = sqrt(max(sigma_Th_100Ppz^2 - T100*(2-T100)*sigma_Th_PPZ^2, 0)) / abs(F_Th_PPZ);
uncertT100_expected = sqrt(45) / 100;   % exact hand derivation: sqrt((136-91)/100^2)
nFail = nFail + checkVal('T100 uncertainty', uncertT100, uncertT100_expected, tol);

%%  verify ratio cancellation: POC-variance formula gives wrong (zero) answer
%   Applying the nested-covariance formula naively with POC variances yields a negative
%   argument because Var_POC is inflated by sigma_R^2; the max guard clamps it to zero.
%   This is the bug the Th-variance formula corrects.
sigma_POC_PPZ    = sqrt(Var_POC_PPZ);
sigma_POC_100Ppz = sqrt(Var_POC_100Ppz);
arg_POC = sigma_POC_100Ppz^2 - T100*(2-T100)*sigma_POC_PPZ^2;
if arg_POC < 0
    fprintf('  pass  POC-variance formula gives negative argument (max guard needed)\n');
else
    fprintf('  FAIL  expected negative argument from POC variances; got %.6g\n', arg_POC);
    nFail = nFail + 1;
end

%%  verify T100 uncertainty > 0 and < 1 (sanity bounds)
if uncertT100 > 0 && uncertT100 < 1
    fprintf('  pass  uncertT100 in (0, 1)\n');
else
    fprintf('  FAIL  uncertT100 = %.6g out of expected range (0, 1)\n', uncertT100);
    nFail = nFail + 1;
end

%%  report
if nFail == 0
    fprintf('testUncertMetrics: ALL PASS\n');
else
    fprintf('testUncertMetrics: %d FAILURE(S)\n', nFail);
end

%%  helper (same contract as testCalcError)
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
