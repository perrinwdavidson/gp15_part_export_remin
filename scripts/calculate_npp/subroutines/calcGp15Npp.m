%%  loop through all measurement locations and interpolate (THIS TAKES HOURS)
%   make sample array ::
Xs = rmmissing([X(:), ...  % [deg]
		Y(:), ...  % [deg] 
	        T(:), ...  % [d]
	        V(:)]);  % [m s-1]
Ts = unique(T);
Nsamp = size(Xs, 1);

%   loop through all depth levels and interpolate (2d) ::
nppStat = []; 
for i = 1 : 1 : NUMSTAT 

	%   get query point ::
	sn = gp15_stations.stationNo(i); 
	xq0 = gp15_stations.longitude(i);
	yq = gp15_stations.latitude(i);
	tq = datenum(gp15_stations.date(i)); 
	tq0 = tq - deltaT; 
	xq = [xq0, yq]; 

	%   get indices of previous time ::
	[~, it1] = min(abs(tq - Ts)); 
	[~, it0] = min(abs(tq0 - Ts)); 

	%   choose depth level ::
	iNpp = []; 
	for j = it0 : 1 : it1
		
		%   get data ::
		Xsj = Xs(T==Ts(j), :);

		%   get sample data distances ::
		deltaX = pdist2(xq, Xsj(:, 1:2));
		idxSamp = (deltaX <= lx);
		x = Xsj(idxSamp, 1:2);
		v = Xsj(idxSamp, 4);
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

		%   store data ::
		iNpp = [iNpp; VhatEstMean];

	end

	%   average and store ::
	nppStat = [nppStat; [sn, tq, mean(iNpp)]];

	%   display ::
	% disp(['Done with station ' num2str(sn)]); 

end

%   make table ::
gp15_npp = array2table(nppStat); 
gp15_npp.Properties.VariableNames = {'Station', 'Sampling_date', 'mmolC'}; 

%   plot ::
figure; 
plot(gp15_stations.latitude, gp15_npp.mmolC, '-ko', 'lineWidth', 2, 'markerFaceColor', 'k'); 
xlabel('\textbf{Latitude [deg N.]}', 'interpreter', 'latex');
ylabel('\textbf{NPP [mmol C m$^{-2}$ d$^{-1}$]}', 'interpreter', 'latex');
set(gca, 'tickLabelInterpreter', 'latex', 'fontSize', 16, 'fontWeight', 'bold', 'lineWidth', 1);
title(['\textbf{CBPM NPP at GP15 Sampling Coordinates}'], 'interpreter', 'latex', 'fontSize', 20);
set(gcf, 'position', [0, 0, 1000, 500]); 
exportgraphics(gcf, [plot_output_basepath 'calcNpp/interp/gp15_npp.png'], 'resolution', 300);

%%  make mean over cruise dates ::
%   get bounds ::
tMin = min(datenum(gp15_stations.date));
tMax = max(datenum(gp15_stations.date));

%   get sample bounds ::
[~, it1] = min(abs(tMax - Ts)); 
[~, it0] = min(abs(tMin - Ts)); 

%   make data ::
npp_mean = mean(NPP(:, :, it0:it1), 3, 'omitnan'); 

%%  end subroutine
