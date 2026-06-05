%% regressionRegionalRatios -- regressing poc:th234 ratios along the gp15 transect
%% ------------------------------------------------------------------------
%% configure
setupModel; 

%% read all data
readRegressionData; 

%% regress data
calcRegionalRegressionPocLsf; 
calcRegionalRegressionPocSsf;  % not needed but good to have
calcRegionalRegressionPnLsf; 
calcRegionalRegressionPnSsf;  % not needed but good to have

%% make ratio model input
makeGp15RegressRatio;  % assume only LSF contributes to flux 

%%  end subroutine
