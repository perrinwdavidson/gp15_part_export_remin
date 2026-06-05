%%  clear and replace
%   setup ::
setupModel; 

%   load ::
load([sim_output_basepath 'gp15Model/modelOutput/gp15_flux.mat'], 'gp15_flux');
load([sim_output_basepath 'gp15Model/modelOutput/gp15_flux100.mat'], 'gp15_flux100');
load([sim_output_basepath 'gp15Model/modelOutput/gp15_fluxPpz.mat'], 'gp15_fluxPpz');
load([sim_output_basepath 'gp15Model/modelOutput/gp15_flux100Ppz.mat'], 'gp15_flux100Ppz');

%%  end subroutine
