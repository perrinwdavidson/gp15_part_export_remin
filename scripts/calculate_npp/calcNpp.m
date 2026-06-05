%%  calcNpp -- calculating npp across the transect
%%  -------------------------------------------------------
%%  configure
setupModel; 

%%  calculate npp
readNppData;
% C1: configureNpp set no variables; deleted
% configureNpp;
calcGp15Npp;
writeNpp; 

%%  write out ::
disp('Done calculating NPP.'); 

%% end routine
