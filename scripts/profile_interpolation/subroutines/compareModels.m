%%  compare models
%   load data ::
w = load([sim_output_basepath 'interpData/upwelling/w_ecco.mat'], 'gp15_w');
i = 1;
clear('gp15_w');
gp15_w{i} = w.gp15_w; 
for dataProduct = {'cglo', 'foam', 'glor', 'oras', 'grep'}
	i = i + 1; 
	w = load([sim_output_basepath 'interpData/upwelling/w_' dataProduct{1} 'Ave.mat'], 'gp15_wAve');
	gp15_w{i} = w.gp15_wAve;  
end

%   get ppz and 100 ::
numModels = length(gp15_w); 
y100mean = zeros(NUMSTAT, 1); 
yPpzmean = zeros(NUMSTAT, 1); 
y100meanGrep = zeros(NUMSTAT, 1); 
yPpzmeanGrep = zeros(NUMSTAT, 1); 
for i = 1 : 1 : numModels
	
	% get data ::
	w100 = gp15_w{i}{gp15_w{i}{:, 4} == 100, :}; 
	[~, idx100, ~] = unique(w100(:, 3));
	[~, idxPpz0] = intersect(gp15_w{i}{:, 4}, gp15_stations.depthPpz); 
	wPpz = table2array(gp15_w{i}(idxPpz0, :)); 
	[~, idxPpz, ~] = unique(wPpz(:, 3));
	
	% make data arrays ::
	x100{i} = w100(idx100, 3);
	y100{i} = w100(idx100, 6) * SEC2DAY;
	n100{i} = sqrt(w100(idx100, 7) ./ w100(:, 8)) * SEC2DAY;
	xPpz{i} = wPpz(idxPpz, 3);
	yPpz{i} = wPpz(idxPpz, 6) * SEC2DAY;
	nPpz{i} = sqrt(wPpz(idxPpz, 7) ./ wPpz(:, 8)) * SEC2DAY;

	% calculate mean ::
	y100mean = y100mean + y100{i}; 
	yPpzmean = yPpzmean + yPpz{i}; 

	% calculate sample mean (sanity check) ::
	if i ~= 1
		y100meanGrep = y100meanGrep + y100{i}; 
		yPpzmeanGrep = yPpzmeanGrep + yPpz{i}; 
	end

end
y100mean = y100mean / numModels; 
yPpzmean = yPpzmean / numModels; 
y100meanGrep = y100meanGrep / (numModels - 1);  % n.b.: grep is mean of other 4 
yPpzmeanGrep = yPpzmeanGrep / (numModels - 1); 

%   calculate l1 norm ::
totDiscrep = []; 
for i = 1 : 1 : (numModels - 1)
	meanDiscrep100{i} = abs(y100{i} - y100mean); 
	meanDiscrepPpz{i} = abs(yPpz{i} - yPpzmean); 
	totDiscrep = [totDiscrep, sum(meanDiscrep100{i}) + sum(meanDiscrepPpz{i})]; 
end

%   plot 
%%% start figure ::
figure; 
tl = tiledlayout(3, 2, 'tileSpacing', 'compact'); 

%%% 100 m ::
nexttile(); 
hold('on'); 
for i = 1 : 1 : length(y100)
	plot(x100{i}, y100{i}, '-o'); 
end
plot(x100{i}, y100mean, '-ok', 'markerFaceColor', 'k');
plot(x100{i}, y100meanGrep, '-ob', 'markerFaceColor', 'b');
xlabel('\textbf{Latitude [deg N.]}', 'interpreter', 'latex');
ylabel('$w_i$ \textbf{[m d$^{-1}$]}', 'interpreter', 'latex');
set(gca, 'tickLabelInterpreter', 'latex', 'fontSize', 16, 'fontWeight', 'bold', 'lineWidth', 1, 'box', 'on');
title('\textbf{PMT Upwelling Velocity at 100 [m]}', 'interpreter', 'latex', 'fontSize', 20);

%%% ppz ::
hold('off'); 
nexttile(); 
hold('on'); 
for i = 1 : 1 : length(yPpz)
	plot(xPpz{i}, yPpz{i}, '-o'); 
end
plot(xPpz{i}, yPpzmean, '-ok', 'markerFaceColor', 'k');
plot(xPpz{i}, yPpzmeanGrep, '-ob', 'markerFaceColor', 'b');
legend('ECCO', 'CGLO', 'FOAM', 'GLOR', 'ORAS', 'GREP', 'Mean w/ ECCO', 'Mean w/o ECCO', 'location', 'eastOutside', 'interpreter', 'latex'); 
hold('off'); 
xlabel('\textbf{Latitude [deg N.]}', 'interpreter', 'latex');
% ylabel('$w_i$ \textbf{[m d$^{-1}$]}', 'interpreter', 'latex');
set(gca, 'tickLabelInterpreter', 'latex', 'fontSize', 16, 'fontWeight', 'bold', 'lineWidth', 1, 'box', 'on');
title('\textbf{PMT Upwelling Velocity at the PPZ}', 'interpreter', 'latex', 'fontSize', 20);

%%% discrepancy at 100 :
nexttile(); 
hold('on'); 
for i = 1 : 1 : length(y100)
	plot(x100{i}, meanDiscrep100{i}, '-o'); 
end
hold('off'); 
xlabel('\textbf{Latitude [deg N.]}', 'interpreter', 'latex');
ylabel('$l^1(w_i)$ \textbf{[m d$^{-1}$]}', 'interpreter', 'latex');
set(gca, 'tickLabelInterpreter', 'latex', 'fontSize', 16, 'fontWeight', 'bold', 'lineWidth', 1, 'box', 'on');
title('\textbf{PMT Absolute Deviation Upwelling Velocity at 100 [m]}', 'interpreter', 'latex', 'fontSize', 20);

%%% discrepancy at the ppz ::
nexttile(); 
hold('on'); 
for i = 1 : 1 : length(yPpz)
	plot(xPpz{i}, meanDiscrep100{i}, '-o'); 
end
hold('off'); 
xlabel('\textbf{Latitude [deg N.]}', 'interpreter', 'latex');
% ylabel('$l^1(w_i)$ \textbf{[m d$^{-1}$]}', 'interpreter', 'latex');
set(gca, 'tickLabelInterpreter', 'latex', 'fontSize', 16, 'fontWeight', 'bold', 'lineWidth', 1, 'box', 'on');
title('\textbf{PMT Absolute Deviation Upwelling Velocity at the PPZ}', 'interpreter', 'latex', 'fontSize', 20);
legend('ECCO', 'CGLO', 'FOAM', 'GLOR', 'ORAS', 'location', 'eastOutside', 'interpreter', 'latex'); 

%%% total discrepancy ::
nexttile(); 
bar(categorical({'ECCO', 'CGLO', 'FOAM', 'GLOR', 'ORAS'}), totDiscrep);
xlabel('\textbf{Model Product}', 'interpreter', 'latex');
ylabel('$\sum_i l^1(w_i)$ \textbf{[m d$^{-1}$]}', 'interpreter', 'latex');
set(gca, 'tickLabelInterpreter', 'latex', 'fontSize', 16, 'fontWeight', 'bold', 'lineWidth', 1, 'box', 'on');
title('\textbf{PMT Total Absolute Deviation Upwelling Velocity}', 'interpreter', 'latex', 'fontSize', 20);

set(gcf, 'position', [0, 0, 1500, 1000]); 

exportgraphics(gcf, [plot_output_basepath 'interpData/upwelling/all_model_transect.png'], 'resolution', 300);

%%  end subroutine
