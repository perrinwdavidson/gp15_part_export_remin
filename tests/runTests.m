%%  runTests - master runner for all automated tests
%   run this script from the project root (or let setupModel set the path).
%   each test prints pass/fail per check, then a summary line.
%   testGoldenFlux requires the full pipeline to have been run at least once. ::

fprintf('\n=== GP15 automated tests ===\n\n');

fprintf('--- testCalcError ---\n');
testCalcError;

fprintf('\n--- testFitPiecewiseRatio ---\n');
testFitPiecewiseRatio;

fprintf('\n--- testUncertMetrics ---\n');
testUncertMetrics;

fprintf('\n--- testGoldenFlux ---\n');
testGoldenFlux;

fprintf('\n=== done ===\n');
