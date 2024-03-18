%% gp15Model - code for kenyon et al. (2024)
%%-------------------------------------------------------------------------
%% configure
clc;
setupModel;

%%  print start
printStartModel;

%%  prepare data
readData;
calcPpz; 
doQc;  
calcNpp;  
calcMld;  
calcEqp;  
calcUpwell; 
interpData;  
calcVertGrad;  
delineateRegions;  
regressRegionalRatio; 
collateData; 

%%  run model
runModel;  

%%  print end
printEndModel;

%% end program
