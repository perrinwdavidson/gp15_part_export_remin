%% make ratios for model input
%  preallocate ::
regressPocRatio = []; 
regressPnRatio = [];
gp15Regions = []; 

%  loop ::
for i = 1 : 1 : size(regions, 1)

	% get stations in region ::
	stats = gp15_stations.stationNo((gp15_stations.stationNo >= regions.first(i)) & (gp15_stations.stationNo <= regions.last(i)));

	% loop ::
	for j = 1 : 1 : length(stats)

		% check in region ::
		sn = stats(j);

		% get data ::
		depthsj = gp15_obsFull.depth(gp15_obsFull.stationNo == sn); 

		% loop through all depths
		for k = 1 : 1 : length(depthsj)

			% get depth ::
			dk = depthsj(k); 

			% get depths ::
			mldi = ratioPocLsf.mldMean(i);
			zstari = ratioPocLsf.zstar(i); 

			% [0, mldMean] ::
			if (dk >= 0) & (dk <= mldi)
				regressPocRatio = [regressPocRatio; [(ratioPocLsf.ratioBeta(i) * mldi) + ratioPocLsf.ratioAlpha(i), ...
							             sqrt(((ratioPocLsf.betaVar(i) * (mldi ^ 2)) + ratioPocLsf.alphaVar(i) + (2 * mldi * ratioPocLsf.coVar(i))) / ratioPocLsf.n(i))]]; 
				regressPnRatio = [regressPnRatio; [(ratioPnLsf.ratioBeta(i) * mldi) + ratioPnLsf.ratioAlpha(i), ...
							           sqrt(((ratioPnLsf.betaVar(i) * (mldi ^ 2)) + ratioPnLsf.alphaVar(i) + (2 * mldi * ratioPnLsf.coVar(i))) / ratioPnLsf.n(i))]]; 

			% (mldMean, zstar] ::
			elseif (dk > ratioPocLsf.mldMean(i)) & (dk <= zstari)
				regressPocRatio = [regressPocRatio; [(ratioPocLsf.ratioBeta(i) * dk) + ratioPocLsf.ratioAlpha(i), ...
							             sqrt(((ratioPocLsf.betaVar(i) * (dk .^ 2)) + ratioPocLsf.alphaVar(i) + (2 * dk * ratioPocLsf.coVar(i))) / ratioPocLsf.n(i))]]; 
				regressPnRatio = [regressPnRatio; [(ratioPnLsf.ratioBeta(i) * dk) + ratioPnLsf.ratioAlpha(i), ...
							           sqrt(((ratioPnLsf.betaVar(i) * (dk .^ 2)) + ratioPnLsf.alphaVar(i) + (2 * dk * ratioPnLsf.coVar(i))) / ratioPnLsf.n(i))]]; 

			% (zstar, bottom] ::
			elseif (dk > zstari) 
				regressPocRatio = [regressPocRatio; [ratioPocLsf.ratioSubsurface(i), ratioPocLsf.ratioSubsurfaceStdErr(i)]];
				regressPnRatio = [regressPnRatio; [ratioPnLsf.ratioSubsurface(i), ratioPnLsf.ratioSubsurfaceStdErr(i)]];
			end

			% make region ::
			gp15Regions = [gp15Regions; [sn, dk, i]]; 

		end

	end

end

%% quality control arrays
regressPocRatio(isnan(gp15_obsFull.pocLarge), :) = NaN; 
regressPnRatio(isnan(gp15_obsFull.pnLarge), :) = NaN; 

%% make final arrays
%  regression ::
regressPocRatio = array2table([gp15_obsFull.stationNo, gp15_obsFull.depth, regressPocRatio]); 
regressPnRatio = array2table([gp15_obsFull.stationNo, gp15_obsFull.depth, regressPnRatio]); 
regressPocRatio.Properties.VariableNames = {'stationNo', 'depth', 'pocTh234RatioRegress', 'uncertPocTh234RatioRegress'};
regressPnRatio.Properties.VariableNames = {'stationNo', 'depth', 'pnTh234RatioRegress', 'uncertPnTh234RatioRegress'};

%  regions ::
gp15Regions = array2table(gp15Regions); 
gp15Regions.Properties.VariableNames = {'stationNo', 'depth', 'region'}; 

%% save
%  .mat ::
save([sim_output_basepath 'regressRegionalRatio/regressPocRatio.mat'], 'regressPocRatio'); 
save([sim_output_basepath 'regressRegionalRatio/regressPnRatio.mat'], 'regressPnRatio'); 
save([sim_output_basepath 'regressRegionalRatio/gp15Regions.mat'], 'gp15Regions'); 

%  .xlsx ::
writetable(regressPocRatio, [sim_output_basepath 'regressRegionalRatio/regressPocRatio.xlsx']); 
writetable(regressPnRatio, [sim_output_basepath 'regressRegionalRatio/regressPnRatio.xlsx']); 
writetable(gp15Regions, [sim_output_basepath 'regressRegionalRatio/gp15Regions.xlsx']); 

%% end subroutine
