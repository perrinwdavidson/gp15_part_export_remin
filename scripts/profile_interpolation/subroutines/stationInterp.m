%%  initialize 
snStart = max(gp15_obs.serialNo); 

%%  spline interpolation per station
for iStat = 1 : 1 : NUMSTAT
    
    %   find station number ::
    statNo = gp15_stations.stationNo(iStat); 
    
    %   find mld and ppz ::
    statMld = gp15_mld.("MLD JAK")(iStat);
    statPpz = gp15_ppz.ppzDepth(iStat); 
    statEqp = gp15_eqp.eqp(iStat);
    statConst = table2array(gp15_stations(iStat, {'depth100', 'depth400'}));
    xq = sort([statMld, statPpz, statPpz + 100, statEqp, statConst]');
    
    %   get station data ::
    statData = gp15_obs(gp15_obs.stationNo == statNo, :);
    statInfo = repmat(statData(1, 1:8), length(xq), 1); 
    statInfo.depth = xq;
    statInfo.serialNo = [snStart + 1 : 1 : snStart + length(xq)]';
    snStart = snStart + length(xq); 

    %	interpolate data ::
    numData = size(statData(:, 9:end), 2);
    statInterpData = NaN(0, numData); 
    for iData = 9 : 2 : size(statData, 2)

	    %    get sample data ::
    	    X = rmmissing([statData.depth, table2array(statData(:, iData)), table2array(statData(:, iData+1))]); 

	    %    get names ::
	    varName = statData.Properties.VariableNames(iData); 
	    varStdName = statData.Properties.VariableNames(iData+1); 

	    %    treat empty ::
	    if isempty(X)

		    %    make update array ::
		    eval(['statInterpData.' varName{1} ' = NaN(length(xq), 1);']); 
		    eval(['statInterpData.' varStdName{1} ' = NaN(length(xq), 1);']); 

	    else

		    %    interpolate ::
		    [statDataInterp, statDataInterpStd] = interp1dError(X(:, 1), X(:, 2), X(:, 3), xq); 	

		    %    make update array ::
		    eval(['statInterpData.' varName{1} ' = statDataInterp;']); 
		    eval(['statInterpData.' varStdName{1} ' = statDataInterpStd;']); 

	    end

    end     

    %   collate data ::
    statDataNew = [statInfo, struct2table(statInterpData)];
    statData = rmmissing(sortrows([statData; statDataNew], 'depth'), 'minNumMissing', numData);

    %   save data ::
    if iStat == 1
	    gp15_obsNew = statData;
    else
	    gp15_obsNew = [gp15_obsNew; statData];
    end

end

%%  average longitude, latitude, and time for all individual stations
for i = 1 : 1 : NUMSTAT

	%   get station number ::
	sn = gp15_stations.stationNo(i); 

	%   replace ::
	gp15_obsNew.longitude(gp15_obsNew.stationNo == sn) = gp15_stations.longitude(i); 
	gp15_obsNew.latitude(gp15_obsNew.stationNo == sn) = gp15_stations.latitude(i); 
	gp15_obsNew.date(gp15_obsNew.stationNo == sn) = gp15_stations.date(i); 

end

%%  write out and clean up
gp15_obs = gp15_obsNew;
clear('gp15_obsNew'); 

%%  save data
%   excel ::
writetable(gp15_obs, [pro_output_basepath 'interpData/observations/gp15_obs.xlsx'])

%   mat ::
save([pro_output_basepath 'interpData/observations/gp15_obs.mat'], 'gp15_obs');

%%  end subroutine
