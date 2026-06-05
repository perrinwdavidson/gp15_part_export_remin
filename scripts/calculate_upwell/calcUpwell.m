%%  calcUpwell - calculating upwelling in a suite of MERCATOR models
%--------------------------------------------------------------------------
%%  set-up environment
setupModel;

%%  load data
loadUpwellData;

%%  configure
configureUpwell;

%%  loop through all model and calculate upwelling
for dataProduct = dataProducts

	%   get product name ::
	dataProduct = dataProduct{1};

	%   read data ::
	readUpwellData;

	%   calculate upwelling ::
	calculateUpwelling;

	%   spatially and temporally average ::
	spatialTimeAverage;

	%   sensitivity: temporal averaging window ::
	sensitivityDeltaDay;

	%   sensitivity: spatial averaging window ::
	sensitivitySpatialAve;

	%   plot ::
	plotUpwell;

	%   save ::
	saveUpwell;

end

%%  process ECCO: load, spatially and temporally average, plot, save
%   ECCO bypasses calculateUpwelling because it provides wo directly rather
%   than u/v fields. it receives the same spatial + temporal averaging
%   pipeline as the MERCATOR products so that kriging training data is
%   consistent in magnitude and smoothing across all six products — the
%   Th-234 35-d integration window applies equally regardless of how w was
%   derived.
dataProduct = 'ecco';

%   load raw ECCO vertical velocity ::
load([pro_output_basepath 'readData/ecco/ecco_wo.mat'], 'ecco_wo');
w = ecco_wo.wo;

%   populate u_mercator interface so that spatialTimeAverage, plotUpwell,
%   sensitivitySpatialAve, and saveUpwell operate without modification ::
u_mercator.latitude  = ecco_wo.latitude;
u_mercator.longitude = ecco_wo.longitude;
u_mercator.depth     = ecco_wo.depth;
u_mercator.time      = ecco_wo.time;
NUMLAT   = length(ecco_wo.latitude);
NUMLON   = length(ecco_wo.longitude);
NUMDEPTH = length(ecco_wo.depth);
NUMTIME  = length(ecco_wo.time);

%   ECCO-specific averaging parameters.
%   spaceResolution = 1.0 deg (ECCO native 1-deg lat-lon grid confirmed from
%   ecco_wo_quart*.cdf: lat/lon spacing = 1.0 deg).
%   timeResolution  = 10  d  (ECCO Wave field is a 10-day average per the
%   ECCO readme; 240 time-step interval at 1-h model step = 10 d, confirmed
%   from ecco_wo_quart*.cdf time diffs). timeAve = floor(35/10) = 3, giving
%   a 30-day moving average — the nearest feasible approximation to the 35-d
%   Th-234 integration window at 10-day resolution. ::
spaceResolution = 1.0;   % deg
% timeResolution  = 30.0;  % d — incorrect: ECCO Wave is 10-day, not monthly
timeResolution  = 10.0;  % d

%   spatially and temporally average ::
spatialTimeAverage;

%   sensitivity: spatial averaging window ::
sensitivitySpatialAve;

%   plot ::
plotUpwell;

%   save ::
saveUpwell;

clear('ecco_wo');

%%  compare ensemble members, their mean, and GREP
compareUpwellModels;

%%  print end
disp('Done calculating upwelling.')

%%  end program
