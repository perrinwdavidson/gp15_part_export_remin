%%  calcPpz -- to calculate ppz depths along the gp15 transect
%%  -------------------------------------------------------
%%  configure
setupModel; 

%%  calculate 
loadPpzData;
calcPpzDepths; 
writePpz; 

%%  print out ::
disp('Done calculating PPZ.');

%%  end routine
