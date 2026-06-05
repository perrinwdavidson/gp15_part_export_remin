%%  calcMld - calculating gp15 MLDs
%--------------------------------------------------------------------------
%%  set-up environment
setupModel;
setModelCoefficients;

%%  configure environment
configureCalcMld;

%%  load data
loadMldData;

%%  mld calculations and plotting
calcMlDepths;

%%  sensitivity: spline smoothing parameter
sensitivitySplineDensity;

%%  export mld data
saveMld;

%% plot mld data 
plotMld;

%% print out 
disp('Done calculating MLDs.')

%% end routine
