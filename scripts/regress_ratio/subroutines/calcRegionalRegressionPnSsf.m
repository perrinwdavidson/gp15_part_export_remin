%% regress data
%  initialize output ::
regress_coeffs = NaN(size(regions, 1), 6);  % of form [beta, alpha, betaVar, alphaVar, coVar, n]
deep_mean = NaN(size(regions, 1), 3);  % of form [mean, stderr, std] 
y_p_ppz = cell(size(regions, 1), 1);
y_p_400 = cell(size(regions, 1), 1);
data_store_ppz = NaN(0, 4); 
data_store_400 = NaN(0, 4); 
sd_ppz = NaN(0, 1); 
sd_400 = NaN(0, 1); 
zstar_region = zeros(size(regions, 1), 1);  
idx_store_ppz = [];
idx_store_400 = [];

%  set bounds ::
BTMDEPTH = 400; 
maxPnSsfPpz = 100;  %3.0;
maxPnSsf400 = 100;  %3.0; 
minPnSsfPpz = 0; %0.05;
minPnSsf400 = 0; %0.05; 

%  calculate average ppz ::
ez_region = zeros(size(regions, 1), 1);  
mld_region = zeros(size(regions, 1), 1);  
for istat = 1 : 1 : size(regions, 1)

	% get bounds ::
	stat_bounds = table2array(regions(istat, 2:3)); 

	% get all data within bounds ::
	idx_stat = find((gp15_stations.stationNo >= stat_bounds(1)) & (gp15_stations.stationNo <= stat_bounds(2))); 
	idx_stat_mld = find((gp15_mld.('Station No') >= stat_bounds(1)) & (gp15_mld.('Station No') <= stat_bounds(2))); 
	
	% get ez data ::
	ez_region(istat) = mean(gp15_stations.depthPpz(idx_stat));
	mld_region(istat) = mean(gp15_mld.('MLD JAK')(idx_stat_mld));

end

%  loop through all stations ::
figure; 
tl = tiledlayout(1, size(regions, 1), 'tileSpacing', 'compact'); 
for istat = 1 : 1 : size(regions, 1)

	% get stations ::
	stat_bounds = table2array(regions(istat, 2:3)); 

	% get all data within bounds ::
	idx_stat = find((gp15_stations.stationNo >= stat_bounds(1)) & (gp15_stations.stationNo <= stat_bounds(2))); 
	data_ppz = NaN(0, 3); 
	data_400 = NaN(0, 3); 
	stat_mld_plot = BTMDEPTH;
	stat_ppz_plot = 0; 
	for idat = 1 : 1 : length(idx_stat)
		stat = idx_stat(idat); 
		stat_mld = gp15_mld.('MLD JAK')(gp15_mld.('Station No') == stat); 
		stat_ppz = gp15_stations.depthPpz(gp15_stations.stationNo == stat);
		if isempty(stat_mld) || isempty(stat_ppz)
			continue
		end
		idx_ppz = find((gp15_obs.stationNo == stat) & (gp15_obs.depth >= stat_mld) & (gp15_obs.depth < stat_ppz) & (~isnan(gp15_obs.pnSmall)) & (gp15_obs.pnTh234RatioSmall <= maxPnSsfPpz) & (gp15_obs.pnTh234RatioSmall >= minPnSsfPpz));
		idx_400 = find((gp15_obs.stationNo == stat) & (gp15_obs.depth >= stat_ppz) & (gp15_obs.depth < BTMDEPTH) & (~isnan(gp15_obs.pnSmall)) & (gp15_obs.pnTh234RatioSmall <= maxPnSsf400) & (gp15_obs.pnTh234RatioSmall >= minPnSsfPpz));
		data_ppz = [data_ppz; [gp15_obs.depth(idx_ppz), gp15_obs.pnTh234RatioSmall(idx_ppz), gp15_obs.uncertPnTh234RatioSmall(idx_ppz)]];
		data_400 = [data_400; [gp15_obs.depth(idx_400), gp15_obs.pnTh234RatioSmall(idx_400), gp15_obs.uncertPnTh234RatioSmall(idx_400)]];
		if stat_mld < stat_mld_plot
			stat_mld_plot = stat_mld; 
		elseif stat_ppz > stat_ppz_plot
			stat_ppz_plot = stat_ppz; 
		end
		idx_store_ppz = [idx_store_ppz; idx_ppz]; 
		idx_store_400 = [idx_store_400; idx_400]; 
	end

	% qc data ::
	data_ppz = rmmissing(data_ppz); 
	data_400 = rmmissing(data_400); 
	data_ppz = sortrows(data_ppz, 1); 
	data_400 = sortrows(data_400, 1); 

	% store ::
	data_store_ppz = [data_store_ppz; [istat*ones(size(data_ppz, 1), 1), data_ppz]];
	data_store_400 = [data_store_400; [istat*ones(size(data_400, 1), 1), data_400]];

	sd_ppz = [sd_ppz; std(data_ppz(:, 2))];
	sd_400 = [sd_400; std(data_400(:, 2))];

	% regress data ::	
	[beta, alpha, q, ci] = ols(data_ppz(:, 1), data_ppz(:, 2)); 
	regress_coeffs(istat, 1:6) = [beta, alpha, ci];
	deep_mean(istat, 1:3) = [mean(data_400(:, 2)), std(data_400(:, 2)) / sqrt(length(data_400(:, 2))), std(data_400(:, 2))];

	% calculate intersection point ::
	if istat == 2
		zstar_region(istat) = ez_region(istat); 
	else
		zstar_region(istat) = (deep_mean(istat, 1) - alpha) / beta;
	end

	% make plotting arrays ::
	depth_ppz = linspace(mld_region(istat), zstar_region(istat), 1000);
	depth_400 = linspace(zstar_region(istat), BTMDEPTH, 1000);

	% calculate plot ::
	y_fit = regress_coeffs(istat, 1)*depth_ppz + regress_coeffs(istat, 2);
	q_interp = interp1(data_ppz(:, 1), q(:, 1), depth_ppz, 'linear', 'extrap');
	y_p_ppz{istat} = [y_fit + q_interp; y_fit - q_interp];

	% plot ::
	nexttile();
	hold('on');

	scatter(data_ppz(:, 2), data_ppz(:, 1), 80, 'k', 'filled');
	
	plot(y_fit, depth_ppz, 'r', 'lineWidth', 2);

	%if istat ~= 1
		plot(y_p_ppz{istat}(1, :), depth_ppz, 'r', 'lineStyle', '--'); 
		plot(y_p_ppz{istat}(2, :), depth_ppz, 'r', 'lineStyle', '--'); 
	%end

	yline(mld_region(istat), 'k', 'lineStyle', '--')
	yline(ez_region(istat), 'k', 'lineStyle', '-.')
	yline(100, 'k', 'lineStyle', ':')

	scatter(data_400(:, 2), data_400(:, 1), 80, 'k', 'filled');
	
	x_tilde = depth_400; 
	plot(repmat(deep_mean(istat, 1), 1, 1000), depth_400, 'r', 'lineWidth', 2);
	y_p_400{istat} = [repmat(deep_mean(istat, 1), 1, 1000) - deep_mean(istat, 2); repmat(deep_mean(istat, 1), 1, 1000) + deep_mean(istat, 2)]';  
	plot(y_p_400{istat}(:, 1), depth_400, 'r', 'lineStyle', '--'); 
	plot(y_p_400{istat}(:, 2), depth_400, 'r', 'lineStyle', '--'); 

	hold('on'); 
	title(region_names(istat), 'interpreter', 'latex');
	set(gca, 'box', 'on', 'yDir', 'reverse', 'xLim', [0, 1], 'yLim', [0, 400], 'tickLabelInterpreter', 'latex'); 

end

%  set labels ::
xlabel(tl, 'SSF PN:$^{234}$Th [$\mu$mol dpm$^{-1}$]', 'interpreter', 'latex'); 
ylabel(tl, 'Depth [m]', 'interpreter', 'latex');
set(gcf, 'position', [0, 0, 1200, 400]);
exportgraphics(gcf, [plot_output_basepath 'regressRegionalRatio/regressionPnSsf.png'], 'resolution', 300);

%% save final result
%  make data ::
ratioPnSsf = array2table([regress_coeffs, deep_mean, mld_region, zstar_region, ez_region]); 
ratioPnSsf = [region_names', ratioPnSsf]; 
ratioPnSsf.Properties.VariableNames = {'region', 'ratioBeta', 'ratioAlpha', 'betaVar', 'alphaVar', 'coVar', 'n', 'ratioSubsurface', 'ratioSubsurfaceStdErr', 'ratioSubsurfaceStd', 'mldMean', 'zstar', 'ppzmean'};

%  save ::
save([sim_output_basepath 'regressRegionalRatio/ratioPnSsf.mat'], 'ratioPnSsf'); 
writetable(ratioPnSsf, [sim_output_basepath 'regressRegionalRatio/ratioPnSsf.xlsx']); 

%% end subroutine
