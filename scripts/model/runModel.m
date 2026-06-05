%%  runModel - model of gp15 thorium-234
%%-------------------------------------------------------------------------
%%  setup environment
setupModel;

%%  load data
loadModelData;

%%  set coefficients
setModelCoefficients;

%%  calculate 234th flux 
for iStat = 1 : 1 : NUMSTAT
    
    doModelCalcs; 
    
end

%%  calculate poc fluxes
calcParticulateFlux; 

%%  calculate error
calcError;  

%%  make ppz, 100m, and 100m below the ppz array
makeDepthArrays; 

%%  save data
saveModelData; 

%%  plot data
plotModel;

%%  calculate statistics
calcStats;  

%%  print out
disp('Done calculating Th-234 and POC fluxes along GP15 transect.');

%%  clear and replace
clearEnvironment; 

%% done with program
