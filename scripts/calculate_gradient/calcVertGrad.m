%% calcVertGrad - vertical gradient calculation
%--------------------------------------------------------------------------
%% setup environment
setupModel;
setModelCoefficients;

%% load data
loadGradData;

%%  calculate gradient
doVertGradCalc;

%%  sensitivity: spline smoothing parameter
% sensitivitySplineGrad;

%% save output
saveVertGradCalc;

%% plot output
plotVertGrad;

%%  print out
disp('Done calculating vertical gradient.');

%% end program
