%% doQc - quality controlling gp15 input data
%%-------------------------------------------------------------------------
%% setup environment
setupModel;

%% read data
doPpzQc;
makeStations;
doGp15Qc;

%% print out 
disp('Done quality controlling data.')

%% end routine
