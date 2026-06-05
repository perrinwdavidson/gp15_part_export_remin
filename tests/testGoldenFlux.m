%%  testGoldenFlux - golden-file regression test for model flux outputs
%   on first run (no golden file): saves current gp15_flux100 as the
%   reference and prints a message.
%   on subsequent runs: compares current output against the golden file.
%   columns checked: th234FluxCumul1d, th234FluxCumul2d_ecco, totalPocFluxError, pocFluxCumul2d_ecco.
%   pass criteria: max absolute difference < 1e-10. ::

setupModel;

tol = 1e-10;
nFail = 0;

goldenPath  = [fileparts(mfilename('fullpath')) '/golden/gp15_flux100_golden.mat'];
currentPath = [sim_output_basepath 'gp15Model/modelOutput/gp15_flux100.mat'];

%%  load current output
if ~isfile(currentPath)
    fprintf('testGoldenFlux: SKIP — %s not found; run the full pipeline first.\n', currentPath);
    return
end
load(currentPath, 'gp15_flux100');

%%  first run: save golden
if ~isfile(goldenPath)
    save(goldenPath, 'gp15_flux100');
    fprintf('testGoldenFlux: golden file created at %s\n', goldenPath);
    fprintf('testGoldenFlux: re-run after any model change to detect regressions.\n');
    return
end

%%  compare against golden
load(goldenPath, 'gp15_flux100');
golden = gp15_flux100;
load(currentPath, 'gp15_flux100');
current = gp15_flux100;

cols = {'th234FluxCumul1d', 'th234FluxCumul2d_ecco', 'totalPocFluxError', 'pocFluxCumul2d_ecco'};
for iCol = 1 : 1 : length(cols)
    col = cols{iCol};
    g   = golden.(col);
    c   = current.(col);
    maxDiff = max(abs(c - g), [], 'all', 'omitnan');
    if maxDiff > tol
        fprintf('  FAIL  %s: max |current - golden| = %.2e (tol = %.2e)\n', col, maxDiff, tol);
        nFail = nFail + 1;
    else
        fprintf('  pass  %s (max diff = %.2e)\n', col, maxDiff);
    end
end

%%  report
if nFail == 0
    fprintf('testGoldenFlux: ALL PASS\n');
else
    fprintf('testGoldenFlux: %d FAILURE(S)\n', nFail);
end
