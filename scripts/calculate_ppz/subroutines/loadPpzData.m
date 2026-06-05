%%  read ppz data
%   ctd ::
load([pro_output_basepath 'readData/gp15/gp15_ctd.mat'], 'gp15_ctd');

%   stations and casts ::
load([pro_output_basepath 'readData/gp15/gp15_stat_cast.mat'], 'gp15_stat_cast');

%   obs ::
load([pro_output_basepath 'readData/gp15/gp15_obs.mat'], 'gp15_obs');

%%  end subroutine
