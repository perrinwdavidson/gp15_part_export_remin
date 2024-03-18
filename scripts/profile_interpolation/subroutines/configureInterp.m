%%  configure
%   get pacific ocean in [0, 400] (m) ::
%%% set bounds (10 deg around cruise) ::
pacificLat = [-60, 58]; 
pacificLon = [129, 293];
pacificDepth = [0, 400];

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
