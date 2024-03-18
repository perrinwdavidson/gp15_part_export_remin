%%  calcMld - calculating gp15 MLDs
%--------------------------------------------------------------------------
%%  set-up environment
setupModel;

%%  configure environment
configureCalcMld;

%%  load data
loadMldData;

%%  mld calculations and plotting
calcMlDepths;

%%  export mld data
saveMld;

%% print out 
disp('Done calculating MLDs.')

%% end routine
