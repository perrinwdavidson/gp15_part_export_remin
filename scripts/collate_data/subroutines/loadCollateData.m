%%  load data
%   gp15 observations ::
load([pro_output_basepath 'interpData/observations/gp15_obs.mat'], 'gp15_obs');

%   velocities ::
w = load([sim_output_basepath 'interpData/upwelling/w_ecco.mat'], 'gp15_w');
gp15_w.ecco = w.gp15_w; 
for dataProduct = {'cglo', 'foam', 'glor', 'oras', 'grep'}
	w = load([sim_output_basepath 'interpData/upwelling/w_' dataProduct{1} 'Ave.mat'], 'gp15_wAve');
	eval(['gp15_w.' dataProduct{1} ' = w.gp15_wAve;']);  
end

%   gradient calculations ::
load([sim_output_basepath 'calcVertGrad/gp15_grad.mat'], 'gp15_grad');

%   regression ratio ::
load([sim_output_basepath 'regressRegionalRatio/regressPocRatio.mat'], 'regressPocRatio'); 
load([sim_output_basepath 'regressRegionalRatio/regressPnRatio.mat'], 'regressPnRatio'); 

%   regions ::
load([sim_output_basepath 'regressRegionalRatio/gp15Regions.mat'], 'gp15Regions'); 

%   ppz data ::
load([pro_output_basepath 'doQc/ppz/gp15_ppz.mat'], 'gp15_ppz');

%   mld data ::
load([sim_output_basepath 'calcMld/gp15_mld.mat'], 'gp15_mld');

%   eqp data ::
load([sim_output_basepath 'calcEqp/gp15_eqp.mat'], 'gp15_eqp');

%   station data:
load([pro_output_basepath 'doQc/stations/gp15_stations.mat'], 'gp15_stations', 'NUMSTAT');

%%  end subroutine
