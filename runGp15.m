%% runGp15 - full analysis code for kenyon and davidson et al. (2026)
%%-------------------------------------------------------------------------
%% configure
addpath(genpath(pwd))
close('all')
clear

%% run full model
gp15Model  % run full pipeline to produce new gp15_flux* files
runTests  % run tests to create new golden file. other tests should pass
runTests  % run tests a second time to confirm all four tests pass. if so, freeze output.

%% end program
