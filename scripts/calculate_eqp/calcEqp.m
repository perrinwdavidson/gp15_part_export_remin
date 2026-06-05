%% calcEqp - equilibrium point calculation
%--------------------------------------------------------------------------
%%  set-up environment
setupModel;
setModelCoefficients;

%%  configure
configureCalcEqp;

%%  load data 
loadEqpData;

%%  calculate eqp depths
calcEqpDepths;

%%  plot eqp
plotEqp;

%% print out
disp('Done calculating EQPs.')

%% end routine
