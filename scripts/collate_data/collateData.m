%%  collectData - collecting and collating all model data inouts
%--------------------------------------------------------------------------
%%  set-up environment 
setupModel;

%%  load data
loadCollateData;

%%  collect data per stations                                   
collateStations; 

%%  save data 
saveCollateData;

%%  print
disp('Done collating data.');

%%  end program
