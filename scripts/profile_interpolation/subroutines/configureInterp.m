%%  configure
%   pacific bounding box for kriging training data.
%
%   latitude [-60, 60]:
%     upper bound 60 N ensures all GP15 stations (northernmost ~59 N, near
%     Kodiak Island) fall inside the training domain and are interpolated
%     rather than extrapolated. the southern bound 60 S extends ~42 deg
%     beyond the southernmost station (Tahiti, ~17.5 S), well past what
%     the 500 km horizontal decorrelation scale makes statistically
%     relevant, but retains the full subtropical and subpolar Pacific for
%     hyperparameter estimation.
%
%   longitude [129, 280]:
%     129 E is in the open western Pacific (east of the Philippine Sea) and
%     gives ~80 deg of training data west of the transect at 150 W.
%     280 E (= 80 W) sits just west of the South American Pacific coast at
%     all latitudes, keeping the domain strictly within the Pacific basin.
%     extending further east past the continental barrier would include
%     atlantic ocean grid cells at latitudes where the americas do not
%     block the domain, introducing a stationarity violation in the
%     hyperparameter optimisation.
%
%   depth [0, 400]:
%     matches BTMDEPTH (the flux-model integration ceiling); all query
%     points are in this range. the matern-3/2 covariance at 400 m
%     separation with GpVerticalScaleM = 100 m is < 0.01, so deep data
%     (> 400 m) has effectively zero influence on predictions. including it
%     would inflate the training set ~15x and distort the fitted vertical
%     scale toward abyssal thermohaline dynamics, which differ structurally
%     from wind-driven near-surface upwelling. ::
pacificLat   = [-60, 60];
pacificLon   = [129, 280];
pacificDepth = [0,   400];

%   conversion ::
SEC2DAY = 60 * 60 * 24; 

%   limits ::
BTMDEPTH = pacificDepth(2);  % m
RESTIME = 35;  % d

%   product names ::
dataProducts = {'ecco', 'cglo', 'foam', 'glor', 'oras', 'grep'}; 

%   spatial resolution ::
spaceResolution = 0.25;  % deg

%%  end subroutine
