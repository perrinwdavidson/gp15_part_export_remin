%%  calcNpp -- calculating npp across the transect
%%  -------------------------------------------------------
%%  configure
setupModel; 

%%  calculate npp
readNppData; 
configureNpp;
calcGp15Npp; 
writeNpp; 

%%  write out ::
disp('Done calculating NPP.'); 

%% end routine
