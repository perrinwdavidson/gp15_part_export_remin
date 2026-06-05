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
%   mean Euclidean inter-station spacing in deg; reduces to ~mean(dlat) for
%   this near-meridional transect (dlon ≈ 0 between consecutive stations).
%   sets the spatial smoothing scale to match the resolvable observation
%   spacing (Bindoff & Wunsch 1992, sigma_a priori, eq. 1). ::
degreeAve = mean(sqrt((gp15_stations.latitude(2:end) - gp15_stations.latitude(1:end-1)) .^ 2 + (gp15_stations.longitude(2:end) - gp15_stations.longitude(1:end-1)) .^ 2));

%   how many days do you want to average over?
% deltaDay = 10;  % d — previous value; inconsistent with Th-234 mean lifetime (~35 d)
deltaDay = 35;  % d — Th-234 mean lifetime: t½/ln(2) = 24.10/0.693 = 34.8 d ≈ 35 d

%   conversion ::
SEC2DAY = 60 * 60 * 24; 
KM2M = 1000;
A = 6378.137 * KM2M;  % Earth radius in metres — passed to deg2km to return metres (overrides km convention)

%% end subroutine
