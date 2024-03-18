%% average duplicate values
for iStat = 1 : 1 : NUMSTAT
    
    % get station numbers ::
    sn = gp15_stations.stationNo(iStat);
    
    % get radio data ::
    data = rmmissing(gp15_obs(gp15_obs.stationNo == sn, :));  % here, we prioritize the radionuclide data
    data = sortrows(data, 'depth');    

    % get radio particulate data ::
    dataRadioPart = gp15_part(gp15_part.stationNo == sn, :);
    dataRadioPart = sortrows(dataRadioPart, 'depth');    
    
    % get particulate data ::
    dataPart = gp15_particles(gp15_particles.station_num == sn, :);
    dataPart = sortrows(dataPart, 'depth_m');    
    
    % average total radio values into array ::
    [~, ia, idx] = unique(data.depth, 'stable');
    dataNew = data(ia, :);
    dataNew.th234 = accumarray(idx, data.th234, [], @mean); 
    dataNew.uncertTh234 = accumarray(idx, data.uncertTh234, [], @(x) sqrt(sum(x .^ 2)) / length(x));  % assume to be standard deviation
    dataNew.u238 = accumarray(idx, data.u238 , [], @mean); 
    dataNew.uncertU238 = accumarray(idx, data.uncertU238, [], @(x) sqrt(sum(x .^ 2)) / length(x)); 
    
    % average radio particulate values into array ::
    [~, ia, idx] = unique(dataRadioPart.depth, 'stable');
    dataRadioPartNew = dataRadioPart(ia, :);
    dataRadioPartNew.th234PartSsf = accumarray(idx, dataRadioPart.th234PartSsf, [], @mean); 
    dataRadioPartNew.uncertTh234PartSsf = accumarray(idx, dataRadioPart.uncertTh234PartSsf, [], @(x) sqrt(sum(x .^ 2)) / length(x));  % assume to be standard deviation
    dataRadioPartNew.th234PartLsf = accumarray(idx, dataRadioPart.th234PartLsf, [], @mean); 
    dataRadioPartNew.uncertTh234PartLsf = accumarray(idx, dataRadioPart.uncertTh234PartLsf, [], @(x) sqrt(sum(x .^ 2)) / length(x));  % assume to be standard deviation

    % average particulate values into array ::
    [~, ia, idx] = unique(dataPart.depth_m, 'stable');
    dataPartNew = dataPart(ia, :);
    dataPartNew.POC_LPT_uM = accumarray(idx, dataPart.POC_LPT_uM, [], @mean); 
    dataPartNew.POC_LPT_uM_std = accumarray(idx, dataPart.POC_LPT_uM_std, [], @(x) sqrt(sum(x .^ 2)) / length(x));  % assume to be standard deviation
    dataPartNew.POC_SPT_uM = accumarray(idx, dataPart.POC_SPT_uM, [], @mean); 
    dataPartNew.POC_SPT_uM_std = accumarray(idx, dataPart.POC_SPT_uM_std, [], @(x) sqrt(sum(x .^ 2)) / length(x));  % assume to be standard deviation
    dataPartNew.PN_LPT_uM = accumarray(idx, dataPart.PN_LPT_uM, [], @mean); 
    dataPartNew.PN_LPT_uM_std = accumarray(idx, dataPart.PN_LPT_uM_std, [], @(x) sqrt(sum(x .^ 2)) / length(x));  % assume to be standard deviation
    dataPartNew.PN_SPT_uM = accumarray(idx, dataPart.PN_SPT_uM, [], @mean); 
    dataPartNew.PN_SPT_uM_std = accumarray(idx, dataPart.PN_SPT_uM_std, [], @(x) sqrt(sum(x .^ 2)) / length(x));  % assume to be standard deviation
    
    % make new array ::
    if iStat == 1
        gp15_obsNew = dataNew;
	gp15_partNew = dataRadioPartNew;
	gp15_particlesNew = dataPartNew;
    else
        gp15_obsNew = [gp15_obsNew; dataNew];
        gp15_partNew = [gp15_partNew; dataRadioPartNew];
        gp15_particlesNew = [gp15_particlesNew; dataPartNew];
    end

end

%%  reassign data
gp15_obs = gp15_obsNew;
gp15_part = gp15_partNew;
gp15_particles = gp15_particlesNew;
clear('gp15_obsNew', 'gp15_partNew', 'gp15_particlesNew'); 

%% end subsubroutine
