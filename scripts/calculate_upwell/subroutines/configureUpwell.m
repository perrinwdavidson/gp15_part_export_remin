%%  set data products
%   what data product do you want to draw from?
%       [1] CGLO - C-GLORS05 from CMCC (It)
%       [2] FOAM - GloSea5 from Met Office (UK)
%       [3] GLOR - GLORYS2V4 from Mercator Ocean (Fr)
%       [4] ORAS - ORAS5 from ECMWF
%       [5] GREP - CMEMS Global Ocean Ensemble Reanalysis product
dataProducts = {'cglo', 'foam', 'glor', 'oras', 'grep'}; 

%%  configure
%   what is resolution of these models?
spaceResolution = 0.25;  % deg
timeResolution = 1.0;  % d

%   what is your spatial averaging bounds (in degrees)?
degreeAve = mean(sqrt((gp15_stations.latitude(2:end) - gp15_stations.latitude(1:end-1)) .^ 2 + (gp15_stations.longitude(2:end) - gp15_stations.longitude(1:end-1)) .^ 2));  % bindoff and wunsch, 1992: sigma_apriori, eq (1)

%   how many days do you want to average over?
deltaDay = 10;  % d, average of ECCO, possibly also try: 35, residence time of th234

%   conversion ::
SEC2DAY = 60 * 60 * 24; 
KM2M = 1000;
A = 6378.137 * KM2M;

%% end subroutine
