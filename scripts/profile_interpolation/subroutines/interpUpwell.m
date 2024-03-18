%% load data
if strcmp(dataProduct, 'ecco')
	load([pro_output_basepath 'readData/ecco/ecco_wo.mat'], 'ecco_wo');
	X0 = ecco_wo.longitude;
	Y0 = ecco_wo.latitude; 
	Z0 = ecco_wo.depth; 
	T0 = datenum(ecco_wo.time);
	V0 = ecco_wo.wo; 
	clear('ecco_wo'); 
else
	fileName = [sim_output_basepath 'calcUpwell/w_' dataProduct '.nc'];
	X0 = ncread(fileName, 'longitude');
	Y0 = ncread(fileName, 'latitude');
	Z0 = ncread(fileName, 'depth');
	T0 = ncread(fileName, 'time');
	V0 = ncread(fileName, 'w');
end

%%  data prep
%%% make coordinates ::
[lon, lat, depth, time] = ndgrid(X0, Y0, Z0, T0); 

%%% get data ::
X = lon((X0 >= pacificLon(1)) & (X0 <= pacificLon(2)), ...
        (Y0 >= pacificLat(1)) & (Y0 <= pacificLat(2)), ...
        (Z0 >= pacificDepth(1)) & (Z0 <= pacificDepth(2)), ...
        :); 
Y = lat((X0 >= pacificLon(1)) & (X0 <= pacificLon(2)), ...
        (Y0 >= pacificLat(1)) & (Y0 <= pacificLat(2)), ...
        (Z0 >= pacificDepth(1)) & (Z0 <= pacificDepth(2)), ...
        :); 
Z = depth((X0 >= pacificLon(1)) & (X0 <= pacificLon(2)), ...
          (Y0 >= pacificLat(1)) & (Y0 <= pacificLat(2)), ...
          (Z0 >= pacificDepth(1)) & (Z0 <= pacificDepth(2)), ...
          :); 
T = time((X0 >= pacificLon(1)) & (X0 <= pacificLon(2)), ...
         (Y0 >= pacificLat(1)) & (Y0 <= pacificLat(2)), ...
         (Z0 >= pacificDepth(1)) & (Z0 <= pacificDepth(2)), ...
         :);
V = V0((X0 >= pacificLon(1)) & (X0 <= pacificLon(2)), ...
       (Y0 >= pacificLat(1)) & (Y0 <= pacificLat(2)), ...
       (Z0 >= pacificDepth(1)) & (Z0 <= pacificDepth(2)), ...
       :);

%% possibly detrend
if strcmp(dataProduct, 'ecco')

	%   detrend data ::
	Vtilde = V; 
	for i = 1 : 1 : size(V, 4)
		Vtilde(:, :, :, i) = V(:, :, :, i) - nanmean(V(:, :, :, i), 'all'); 	
	end

	%   make array output ::
	pacificData = rmmissing([X(:), Y(:), T(:), V(:), Vtilde(:)]); 

	%   save for statistical determination of decorrelation scales
	fname = [sim_output_basepath 'interpData/decorrelation/ecco_pacific.nc'];
	delete(fname); 
	vname = 'pacific_wo';
	nccreate(fname, vname, 'dimensions', {'observations', size(pacificData, 1), 'variables', size(pacificData, 2)}, 'fillValue', 'disable');
	ncwrite(fname, vname, pacificData);

end

%%  calculate decorrelation scales by running `calculate_decorrelation.jl'
%   load in ::
lxRaw = readmatrix([sim_output_basepath 'interpData/decorrelation/spatial_decorrelation.csv']);
ltRaw = readmatrix([sim_output_basepath 'interpData/decorrelation/temporal_decorrelation.csv']);

%   appropriate round ::
lx = km2deg(round(lxRaw(2, 1),  1 + int64(floor(log10(lxRaw(2, 2)))), 'significant'));  % [deg]
lt = round(ltRaw(2, 1),  1 + int64(floor(log10(ltRaw(2, 2)))), 'significant');  % [d]

%   set initial length bounds ::
lxCurrent = lx;
ltCurrent = RESTIME; 
lzCurrent = 10; 

%   set tolerances ::
XTOL = 1E-10;
TTOL = 1E-1;
ZTOL = 1E-2;

%%  loop through all measurement locations and interpolate (THIS TAKES HOURS)
%   make sample array ::
Xs = rmmissing([X(:), ...  % [deg]
		Y(:), ...  % [deg] 
		Z(:), ...  % [m] 
	        T(:), ...  % [d]
	        V(:)]);  % [m s-1]
Xs(Xs(:, 5) == 0, :) = [];  % quality control routine
Nsamp = size(Xs, 1);

%   loop through all depth levels and interpolate (2d) ::
for i = 1 : 1 : NUMSTAT 

	%   get query point ::
	sn = gp15_stations.stationNo(i); 
	xq0 = gp15_stations.longitude(i) + 360; 
	yq = gp15_stations.latitude(i);
	xq = [xq0, yq]; 

	%   choose depth level ::
	for j = 1 : 1 : length(Z0(Z0 <= BTMDEPTH))			
		
		%   chose time level ::
		for k = 1 : 1 : length(T0)

			%   get data ::
			Xsjk = squeeze(X(:, :, j, k));
			Ysjk = squeeze(Y(:, :, j, k));
			Vsjk = squeeze(V(:, :, j, k)); 
			XsIn = double(rmmissing([Xsjk(:), Ysjk(:), Vsjk(:)]));	
			XsIn(XsIn(:, 3) == 0, :) = [];  % quality control routine

			%   possibly interpolate ::
			if isempty(XsIn)

				%   make addition ::
				wIntAdd = [sn, xq, double(Z0(j)), datenum(T0(k)), NaN, NaN, NaN, 0];

			else
			
				%   get sample data distances ::
				deltaX = pdist2(xq, XsIn(:, 1:2));
				idxSamp = (deltaX <= lx);
				x = XsIn(idxSamp, 1:2);
				v = XsIn(idxSamp, 3);
				vVar = zeros(size(v)); 
				nx = length(x); 

				%   objectively map ::
				omFlag = true;
				lxPrime = lxCurrent; 
				num = 1;
				while omFlag
					lxPrime = lxPrime / num;
					try
						[VhatEstMean, VhatVar, lxEst] = objectiveMapping(x, v, xq, lxPrime, vVar);
						omFlag = false;

					end
					% disp('Reducing lx.'); 
					num = num + 1; 
					if lxPrime < XTOL
						error('Horizontal length to small.')
					end
				end
				lxCurrent = lxPrime;  % always the same

				%   make array ::
				wIntAdd = [sn, xq, double(Z0(j)), datenum(T0(k)), VhatEstMean, VhatVar, lxEst, nx];

			end

			%   store data ::
			if (i == 1) && (j == 1) && (k == 1)
				wIntLatLon = wIntAdd;
			else
				wIntLatLon = [wIntLatLon; wIntAdd];
			end

		end

	end

	%   display ::
	% disp(['Done with station ' num2str(sn)]); 

end

%%% plot ::
figure; 

tl = tiledlayout(2, 1, 'tileSpacing', 'compact'); 

nexttile(); 
scatter3(wIntLatLon(:, 3), wIntLatLon(:, 5) - min(wIntLatLon(:, 5)), wIntLatLon(:, 4), 50, wIntLatLon(:, 6) * SEC2DAY, 'filled'); 
xlabel('\textbf{Latitude [deg N.]}', 'interpreter', 'latex');
ylabel('\textbf{Time [d]}', 'interpreter', 'latex');
zlabel('\textbf{Depth [m]}', 'interpreter', 'latex');
set(gca, 'zDir', 'reverse', 'tickLabelInterpreter', 'latex', 'fontSize', 16, 'fontWeight', 'bold', 'lineWidth', 1);
title(['\textbf{PMT Upwelling Velocity at GP15 Sampling Coordinates (' upper(dataProduct) ')}'], 'interpreter', 'latex', 'fontSize', 20);
cb = colorbar; 
ylabel(cb, '$w$ [m d$^{-1}$]', 'interpreter', 'latex');
set(cb, 'tickLabelInterpreter', 'latex', 'fontSize', 16, 'fontWeight', 'bold', 'lineWidth', 1); 

nexttile(); 
scatter3(wIntLatLon(:, 3), wIntLatLon(:, 5) - min(wIntLatLon(:, 5)), wIntLatLon(:, 4), 50, sqrt(wIntLatLon(:, 7) ./ wIntLatLon(:, 9)) * SEC2DAY, 'filled'); 
xlabel('\textbf{Latitude [deg N.]}', 'interpreter', 'latex');
ylabel('\textbf{Time [d]}', 'interpreter', 'latex');
zlabel('\textbf{Depth [m]}', 'interpreter', 'latex');
set(gca, 'zDir', 'reverse', 'tickLabelInterpreter', 'latex', 'fontSize', 16, 'fontWeight', 'bold', 'lineWidth', 1);
title(['\textbf{PMT Upwelling Velocity Standard Error at GP15 Sampling Coordinates (' upper(dataProduct) ')}'], 'interpreter', 'latex', 'fontSize', 20);
cb = colorbar; 
ylabel(cb, '$\sigma_{\hat{w}}$ [m d$^{-1}$]', 'interpreter', 'latex'); 
set(cb, 'tickLabelInterpreter', 'latex', 'fontSize', 16, 'fontWeight', 'bold', 'lineWidth', 1); 

set(gcf, 'position', [0, 0, 1000, 1000]); 

exportgraphics(gcf, [plot_output_basepath 'interpData/upwelling/' dataProduct '_sampleLocations.png'], 'resolution', 300);

%%% save data ::
save([sim_output_basepath 'interpData/wIntLatLon_' dataProduct '.mat'], 'wIntLatLon'); 

%   loop through all times and interpolate to time (1d) ::
if strcmp(dataProduct, 'ecco')  

	for i = 1 : 1 : NUMSTAT

		%   get query point ::
		sn = gp15_stations.stationNo(i); 
		xq0 = gp15_stations.longitude(i) + 360; 
		yq = gp15_stations.latitude(i);
		xq = [xq0, yq]; 
		tq = datenum(gp15_stations.date(i));

		%   get sample data ::
		Xs = wIntLatLon(wIntLatLon(:, 1) == sn, :); 

		%   get sample depths ::
		zq = unique(Xs(:, 4)); 

		%   choose depth ::
		for j = 1 : 1 : length(zq)

			%   get depth data ::
			xs = Xs(Xs(:, 4) == zq(j), :); 

			%   get sample data distances ::
			deltaT = pdist2(tq, xs(:, 5));
			idxSamp = (deltaT <= lt);
			x = xs(idxSamp, 5);
			v = xs(idxSamp, 6);
			vVar = xs(idxSamp, 7);
			nx = length(x) + sum(Xs(idxSamp, 9)); 

			%   objectively map ::
			omFlag = true;
			ltPrime = ltCurrent; 
			num = 1;
			while omFlag
				ltPrime = ltPrime / num;
				try
					[VhatEstMean, VhatVar, ltEst] = objectiveMapping(x, v, tq, ltPrime, vVar);
					omFlag = false;

				end
				% disp('Reducing lt.'); 
				num = num + 1; 
				if ltPrime < TTOL
					error('Time length to small.')
				end
			end

			%   make array ::
			wIntAdd = [sn, xq, zq(j), tq, VhatEstMean, VhatVar, ltEst, nx];

			%   store data ::
			if (i == 1) && (j == 1) 
				wIntTime = wIntAdd;
			else
				wIntTime = [wIntTime; wIntAdd];
			end

		end

		%   display ::
		% disp(['Done with station ' num2str(sn)]); 

	end

	%%% plot ::
	figure; 

	tl = tiledlayout(2, 1, 'tileSpacing', 'compact'); 

	nexttile(); 
	scatter3(wIntTime(:, 3), wIntTime(:, 5) - min(wIntTime(:, 5)), wIntTime(:, 4), 200, wIntTime(:, 6) * SEC2DAY, 'filled'); 
	xlabel('\textbf{Latitude [deg N.]}', 'interpreter', 'latex');
	ylabel('\textbf{Time [d]}', 'interpreter', 'latex');
	zlabel('\textbf{Depth [m]}', 'interpreter', 'latex');
	set(gca, 'zDir', 'reverse', 'tickLabelInterpreter', 'latex', 'fontSize', 16, 'fontWeight', 'bold', 'lineWidth', 1);
	title(['\textbf{PMT Upwelling Velocity at GP15 Sampling Dates (' upper(dataProduct) ')}'], 'interpreter', 'latex', 'fontSize', 20);
	cb = colorbar; 
	ylabel(cb, '$w$ [m d$^{-1}$]', 'interpreter', 'latex');
	set(cb, 'tickLabelInterpreter', 'latex', 'fontSize', 16, 'fontWeight', 'bold', 'lineWidth', 1); 

	nexttile(); 
	scatter3(wIntTime(:, 3), wIntTime(:, 5) - min(wIntTime(:, 5)), wIntTime(:, 4), 200, sqrt(wIntTime(:, 7) ./ wIntTime(:, 9)) * SEC2DAY, 'filled'); 
	xlabel('\textbf{Latitude [deg N.]}', 'interpreter', 'latex');
	ylabel('\textbf{Time [d]}', 'interpreter', 'latex');
	zlabel('\textbf{Depth [m]}', 'interpreter', 'latex');
	set(gca, 'zDir', 'reverse', 'tickLabelInterpreter', 'latex', 'fontSize', 16, 'fontWeight', 'bold', 'lineWidth', 1);
	title(['\textbf{PMT Upwelling Velocity Standard Error at GP15 Sampling Dates (' upper(dataProduct) ')}'], 'interpreter', 'latex', 'fontSize', 20);
	cb = colorbar; 
	ylabel(cb, '$\sigma_{\hat{w}}$ [m d$^{-1}$]', 'interpreter', 'latex'); 
	set(cb, 'tickLabelInterpreter', 'latex', 'fontSize', 16, 'fontWeight', 'bold', 'lineWidth', 1); 

	set(gcf, 'position', [0, 0, 1000, 1000]); 

	exportgraphics(gcf, [plot_output_basepath 'interpData/upwelling/' dataProduct '_sampleTime.png'], 'resolution', 300);

	%%% save data ::
	save([sim_output_basepath 'interpData/wIntTime_' dataProduct '.mat'], 'wIntTime'); 

end

%   loop through all stations and interpolate to depth (1d) ::
%%% set averaging time ::
TIMEAVE = 16;  % d, same time as NPP analysis, other options include: RESTIME 

%%% loop through all ::
for i = 1 : 1 : NUMSTAT

	%   get query point ::
	sn = gp15_stations.stationNo(i); 
	xq0 = gp15_stations.longitude(i) + 360; 
	yq = gp15_stations.latitude(i);
	xq = [xq0, yq]; 
	zq = gp15_obs.depth((gp15_obs.stationNo == sn) & (gp15_obs.depth <= BTMDEPTH)); 
	tq = datenum(gp15_stations.date(i));

	%   get sample data ::
	if strcmp(dataProduct, 'ecco')
		Xs = [wIntLatLon(wIntLatLon(:, 1) == sn, :); wIntTime(wIntTime(:, 1) == sn, :);];
	else
		Xs = wIntLatLon(wIntLatLon(:, 1) == sn, :);  % as MERCATOR at daily resolution
	end
	Xs = rmmissing(Xs); 

	%   get unique time array ::
	Ts = unique(Xs((Xs(:, 5) >= tq - TIMEAVE) & (Xs(:, 5) <= tq), 5)); 

	%   loop through all time arrays ::
	for j = 1 : 1 : length(Ts)

		%   get time ::
		tj = Ts(j); 

		%   get sample data distances ::
		idxj = Xs(:, 5) == tj; 
		x = Xs(idxj, 4);
		v = Xs(idxj, 6);
		vVar = Xs(idxj, 7);
		nx = length(x) + sum(Xs(idxj, 9)); 

		%   objectively map ::
		omFlag = true;
		lzPrime = lzCurrent; 
		num = 1;
		while omFlag
			lzPrime = lzPrime / num;
			try
				[VhatEstMean, VhatVar, lzEst] = objectiveMapping(x, v, zq, lzPrime, vVar);
				omFlag = false;

			end
			% disp('Reducing lz.'); 
			num = num + 1; 
			if sum(lzPrime < ZTOL) > 0
				error('Depth length scale to small.')
			end
		end
		lzCurrent = lzPrime; 

		%   store data ::
		if j == 1
			wInt = VhatEstMean;
			wIntVar = VhatVar;
			wIntN = repmat(nx, length(zq), 1); 
		else
			wInt = [wInt, VhatEstMean];
			wIntVar = [wIntVar, VhatVar];
			wIntN = [wIntN, repmat(nx, length(zq), 1)];
		end

	end

	%   calculate mean ::
	wMean = mean(wInt, 2); 
	wVar = sum(wIntVar, 2) / (length(wIntVar) ^ 2); 
	wN = sum(wIntN, 2); 

	%   collate ::
	xqin = repmat([sn, xq], length(zq), 1); 
	tqin = repmat(tq, length(zq), 1);
	gp15_wAdd = [xqin, zq, tqin, wMean, wVar, wN];

	%   store data ::
	if i == 1
		gp15_w = gp15_wAdd;
	else
		gp15_w = [gp15_w; gp15_wAdd];
	end

	%   display ::
	% disp(['Done with station ' num2str(sn)]); 

end

%%  possibly bin data
if ~strcmp(dataProduct, 'ecco')

	%   copy data ::
	gp15_wAve = gp15_w;

	%   make bins ::
	depthBinEdges = [0 : 100 : BTMDEPTH] + 1;

	%   get indices ::
	[~, k] = histc(gp15_w(:, 4), depthBinEdges); 
	ku = unique(k);
	nku = length(ku);

	%   start plot ::
	figure; 
	tl = tiledlayout(floor(sqrt(length(depthBinEdges))), floor(sqrt(length(depthBinEdges))) + mod(nku, 2), 'tileSpacing', 'compact'); 

	%   loop ::
	for i = 1 : 1 : nku
		
		%   get data ::
		ik = ku(i);
		idxk = k==ik;
		s = gp15_wAve(idxk, 1);
		x = gp15_wAve(idxk, 3); 
		z = gp15_wAve(idxk, 4);
		v = gp15_wAve(idxk, 6);
		vvar = gp15_wAve(idxk, 7);
		xn = gp15_wAve(idxk, 8);

		%   get stations ::
		sn = unique(s); 

		%   average by station ::
		for j = 1 : 1 : length(sn)

			% get data ::
			xs = x(s == sn(j)); 
			vs = v(s == sn(j));
			vsvar = vvar(s == sn(j));
			xsn = xn(s == sn(j));

			% calculate average, propagate variance, and count ::
			vsmean = mean(vs); 
			vsmeanvar = sum(vsvar) / (length(vsvar) ^ 2); 
			xsnt = sum(xsn); 

			% store ::
			if j == 1
				xmean = mean(xs); 
				vmean = vsmean; 
				vvarmean = vsmeanvar; 
				nmean = xsnt; 
			else
				xmean = [xmean, mean(xs)]; 
				vmean = [vmean, vsmean]; 
				vvarmean = [vvarmean, vsmeanvar]; 
				nmean = [nmean, xsnt];
			end

			% put into array ::
			gp15_wAve(idxk & (gp15_wAve(:, 1) == sn(j)), 6) = vsmean; 	
			gp15_wAve(idxk & (gp15_wAve(:, 1) == sn(j)), 7) = vsmeanvar; 	
			gp15_wAve(idxk & (gp15_wAve(:, 1) == sn(j)), 8) = xsnt; 	

		end

		%   make data :;
		wAveAdd = [sn, repmat(mean(depthBinEdges(i:i+1)), length(sn), 1), vmean', vvarmean', nmean'];
		
		%   store ::
		if i == 1
			wAve = wAveAdd;
		else
			wAve = [wAve; wAveAdd];
		end

		%   plot
		nexttile(); 
		hold('on');  
		plot(xmean, vmean * SEC2DAY, '--k','lineWidth', 1.5)
		err1 = errorbar(xmean, vmean * SEC2DAY, sqrt(vvarmean ./ nmean), sqrt(vvarmean ./ nmean), '--k', 'lineWidth', 1.5); 
		err1.Color = [0 0 0];                            
		err1.LineStyle = 'none';
		yline(0, ':k', 'lineWidth', 1.5); 
		hold('off'); 
		title(['(' num2str(depthBinEdges(i)-1) ', ' num2str(depthBinEdges(i+1)-1) '] Depth Bin'], 'interpreter', 'latex', 'fontSize', 16); 
		legend('$w$', '', '$\epsilon_{w}$', 'interpreter', 'latex', 'fontSize', 16, 'location', 'southWest'); 
		set(gca, 'box', 'on', 'tickLabelInterpreter', 'latex', 'fontSize', 16, 'fontWeight', 'bold', 'lineWidth', 1);

	end

	%   finish plots ::
	xlabel(tl, '\textbf{Latitude (deg. N)}', 'interpreter', 'latex', 'fontSize', 20);
	ylabel(tl, '\textbf{Upwelling velocity, $w$ (m day$^{-1}$)}', 'interpreter', 'latex', 'fontSize', 20);
	set(gcf, 'units', 'inches', 'position', [1.5, 1.5, 16, 7], 'paperUnits', 'inches', 'paperSize', [16, 7]);
	exportgraphics(gcf, [plot_output_basepath 'interpData/upwelling/' dataProduct '_transect_depths.png'], 'resolution', 300);

else

	%%% get data ::
	w100 = gp15_w(gp15_w(:, 4) == 100, :); 
	[~, idx100, ~] = unique(w100(:, 3));
	[~, idxPpz0] = intersect(gp15_w(:, 4), gp15_stations.depthPpz); 
	wPpz = gp15_w(idxPpz0, :); 
	[~, idxPpz, ~] = unique(wPpz(:, 3));

	%%% make data ::
	x100 = w100(idx100, 3);
	y100 = w100(idx100, 6) * SEC2DAY;
	n100 = sqrt(w100(idx100, 7) ./ w100(:, 8)) * SEC2DAY;
	xPpz = wPpz(idxPpz, 3);
	yPpz = wPpz(idxPpz, 6) * SEC2DAY;
	nPpz = sqrt(wPpz(idxPpz, 7) ./ wPpz(:, 8)) * SEC2DAY;

	%%% plot contour ::
	figure; 

	tl = tiledlayout(2, 1, 'tileSpacing', 'compact'); 

	nexttile(); 
	scatter(gp15_w(:, 3), gp15_w(:, 4), 200, gp15_w(:, 6) * SEC2DAY, 'filled'); 
	xlabel('\textbf{Latitude [deg N.]}', 'interpreter', 'latex');
	ylabel('\textbf{Depth [m]}', 'interpreter', 'latex');
	set(gca, 'yDir', 'reverse', 'tickLabelInterpreter', 'latex', 'fontSize', 16, 'fontWeight', 'bold', 'lineWidth', 1);
	title(['\textbf{PMT Upwelling Velocity (' upper(dataProduct) ')}'], 'interpreter', 'latex', 'fontSize', 20);
	cb = colorbar; 
	ylabel(cb, 'w [m d$^{-1}$]', 'interpreter', 'latex');
	set(cb, 'tickLabelInterpreter', 'latex', 'fontSize', 16, 'fontWeight', 'bold', 'lineWidth', 1); 
	caxis([-5 5]);

	nexttile(); 
	scatter(gp15_w(:, 3), gp15_w(:, 4), 200, sqrt(gp15_w(:, 7) ./ gp15_w(:, 8)) * SEC2DAY, 'filled'); 
	xlabel('\textbf{Latitude [deg N.]}', 'interpreter', 'latex');
	ylabel('\textbf{Depth [m]}', 'interpreter', 'latex');
	set(gca, 'yDir', 'reverse', 'tickLabelInterpreter', 'latex', 'fontSize', 16, 'fontWeight', 'bold', 'lineWidth', 1);
	title(['\textbf{PMT Upwelling Velocity Standard Error (' upper(dataProduct) ')}'], 'interpreter', 'latex', 'fontSize', 20);
	cb = colorbar; 
	ylabel(cb, '$\sigma_{\hat{w}}$ [m d$^{-1}$]', 'interpreter', 'latex'); 
	set(cb, 'tickLabelInterpreter', 'latex', 'fontSize', 16, 'fontWeight', 'bold', 'lineWidth', 1); 
	caxis([-0.5 0.5]);

	set(gcf, 'position', [0, 0, 1000, 1000]); 

	exportgraphics(gcf, [plot_output_basepath 'interpData/upwelling/' dataProduct '_transect.png'], 'resolution', 300);

	%%% plot line plot
	figure; 

	hold('on');  

	plot(x100, y100, '--k','lineWidth', 1.5)
	err1 = errorbar(x100, y100, n100,  n100, '--k', 'lineWidth', 1.5); 
	err1.Color = [0 0 0];                            
	err1.LineStyle = 'none';

	plot(xPpz, yPpz, '-k', 'lineWidth', 1.5)
	err2 = errorbar(xPpz, yPpz, nPpz,  nPpz, '-k', 'lineWidth', 1.5); 
	err2.Color = [0 0 0];                            
	err2.LineStyle = 'none';

	yline(0, ':k', 'lineWidth', 1.5); 

	hold('off'); 

	xlabel('\textbf{Latitude (deg. N)}', 'interpreter', 'latex', 'fontSize', 20);
	ylabel('\textbf{Upwelling velocity, $w$ (m day$^{-1}$)}', 'interpreter', 'latex', 'fontSize', 20);
	legend('$w_{100}$', '', '$w_{PPZ}$', '$\epsilon_{w}$', 'interpreter', 'latex', 'fontSize', 16, 'location', 'eastOutside'); 
	    
	set(gca, 'box', 'on', 'tickLabelInterpreter', 'latex', 'fontSize', 16, 'fontWeight', 'bold', 'lineWidth', 1);

	set(gcf, 'units', 'inches', 'position', [1.5, 1.5, 16, 7], 'paperUnits', 'inches', 'paperSize', [16, 7]);

	exportgraphics(gcf, [plot_output_basepath 'interpData/upwelling/' dataProduct '_transect_depths.png'], 'resolution', 300);

end

%%  save data
%   raw ::
if strcmp(dataProduct, 'ecco')

	% make data ::
	gp15_w = [gp15_w(:, 1:4), datenum(gp15_w(:, 5)), gp15_w(:, 6:7), sqrt(gp15_w(:, 7) ./ gp15_w(:, 8)), gp15_w(:, 8)];
	gp15_w = array2table(gp15_w, 'variableNames', {'stationNo', 'longitude', 'latitude', 'depth', 'time', 'w', 'wVar', 'wErr', 'wN'}); 

	% save data ::
	writetable(gp15_w, [sim_output_basepath 'interpData/upwelling/w_' dataProduct '.xlsx']); 
	save([sim_output_basepath 'interpData/upwelling/w_' dataProduct '.mat'], 'gp15_w'); 

else

	% make data ::
	gp15_wAve = [gp15_wAve(:, 1:4), datenum(gp15_wAve(:, 5)), gp15_wAve(:, 6:7), sqrt(gp15_wAve(:, 7) ./ gp15_wAve(:, 8)), gp15_wAve(:, 8)];
	gp15_wAve = array2table(gp15_wAve, 'variableNames', {'stationNo', 'longitude', 'latitude', 'depth', 'time', 'w', 'wVar', 'wErr', 'wN'}); 

	% save data ::
	writetable(gp15_wAve, [sim_output_basepath 'interpData/upwelling/w_' dataProduct 'Ave.xlsx']); 
	save([sim_output_basepath 'interpData/upwelling/w_' dataProduct 'Ave.mat'], 'gp15_wAve'); 

end

%%  done with subroutine
