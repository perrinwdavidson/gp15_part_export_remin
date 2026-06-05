%%  compareModels - compare kriged upwelling velocity across all six products
%   loads the kriged w for all six products at 100 m and PPZ; produces two
%   figures, each with two stacked panels (100 m, PPZ):
%
%   figure 1 — full ensemble:
%     CMEMS members (CGLO, FOAM, GLOR, ORAS) — gray, distinct line styles
%     CMEMS mean (CGLO+FOAM+GLOR+ORAS)       — thick black
%     GREP                                    — thinner red
%     ECCO                                    — dark blue
%     all-model mean (five independent)       — thick dark green
%
%   figure 2 — FOAM and GREP excluded (sensitivity):
%     CMEMS members (CGLO, GLOR, ORAS)        — gray, same line styles as fig 1
%     CMEMS mean (CGLO+GLOR+ORAS)             — thick black
%     ECCO                                    — dark blue
%     all-model mean (ECCO+CGLO+GLOR+ORAS)    — thick dark green
%--------------------------------------------------------------------------

%%  load kriged upwelling for all six products
w = load([sim_output_basepath 'interpData/upwelling/w_ecco.mat'], 'gp15_w');
i = 1;
clear('gp15_w');
gp15_w{i} = w.gp15_w;
for dataProduct = {'cglo', 'foam', 'glor', 'oras', 'grep'}
	i = i + 1;
	w = load([sim_output_basepath 'interpData/upwelling/w_' dataProduct{1} 'Ave.mat'], 'gp15_wAve');
	gp15_w{i} = w.gp15_wAve;
end

%%  product indices
%   GREP (i=6) is the arithmetic mean of CGLO/FOAM/GLOR/ORAS — not an
%   independent member. ensemble mean uses only IDX_INDEP (n=5). ::
IDX_ECCO  = 1;
IDX_CMEMS = 2 : 5;   % CGLO, FOAM, GLOR, ORAS
IDX_GREP  = 6;
IDX_INDEP = 1 : 5;   % ECCO + four CMEMS members
numModels = length(gp15_w);

%%  extract w at 100 m and PPZ for each product
x100 = cell(numModels, 1);
y100 = cell(numModels, 1);
xPpz = cell(numModels, 1);
yPpz = cell(numModels, 1);
for i = 1 : 1 : numModels

	%   100 m ::
	w100          = gp15_w{i}{gp15_w{i}{:, 4} == 100, :};
	[~, idx100, ~] = unique(w100(:, 3));
	x100{i}       = w100(idx100, 3);
	y100{i}       = w100(idx100, 6) * SEC2DAY;

	%   PPZ ::
	[~, idxPpz0]   = intersect(gp15_w{i}{:, 4}, gp15_stations.depthPpz);
	wPpz_i         = table2array(gp15_w{i}(idxPpz0, :));
	[~, idxPpz, ~] = unique(wPpz_i(:, 3));
	xPpz{i}        = wPpz_i(idxPpz, 3);
	yPpz{i}        = wPpz_i(idxPpz, 6) * SEC2DAY;

end

%%  ensemble means
%   CMEMS mean (CGLO, FOAM, GLOR, ORAS only) ::
y100cmemsMean = mean(cat(2, y100{IDX_CMEMS}), 2, 'omitnan');
yPpzcmemsMean = mean(cat(2, yPpz{IDX_CMEMS}), 2, 'omitnan');

%   all-model mean (five independent products: ECCO + CMEMS) ::
y100allMean = mean(cat(2, y100{IDX_INDEP}), 2, 'omitnan');
yPpzallMean = mean(cat(2, yPpz{IDX_INDEP}), 2, 'omitnan');

%%  styling
cmems_labels     = {'CGLO', 'FOAM', 'GLOR', 'ORAS'};
cmems_lineStyles = {'-', '--', '-.', ':'};
cmems_color      = [0.65, 0.65, 0.65];
ecco_color       = [0.08, 0.17, 0.55];
green_color      = [0.05, 0.45, 0.15];

%%  plot
figure;
tl = tiledlayout(2, 1, 'tileSpacing', 'compact');

%%% upwelling at 100 m ::
nexttile();
hold('on');
for iC = 1 : 1 : length(IDX_CMEMS)
	plot(x100{IDX_CMEMS(iC)}, y100{IDX_CMEMS(iC)}, cmems_lineStyles{iC}, ...
	     'color', cmems_color, 'lineWidth', 1, 'displayName', cmems_labels{iC});
end
plot(x100{IDX_ECCO}, y100cmemsMean,      '-k',                       'lineWidth', 2.5, 'displayName', 'CMEMS Mean');
plot(x100{IDX_GREP}, y100{IDX_GREP},     '-',  'color', 'r',         'lineWidth', 1.5, 'displayName', 'GREP');
plot(x100{IDX_ECCO}, y100{IDX_ECCO},     '-',  'color', ecco_color,  'lineWidth', 2,   'displayName', 'ECCO');
plot(x100{IDX_ECCO}, y100allMean,        '-',  'color', green_color, 'lineWidth', 2.5, 'displayName', 'All-Model Mean');
xline(0, '-k', 'lineWidth', 0.5, 'handleVisibility', 'off');
hold('off');
ylabel('$w$ \textbf{[m d$^{-1}$]}', 'interpreter', 'latex', 'fontSize', 16);
title('\textbf{Kriged Upwelling Velocity at 100 m}', 'interpreter', 'latex', 'fontSize', 18);
set(gca, 'tickLabelInterpreter', 'latex', 'fontSize', 14, 'fontWeight', 'bold', 'lineWidth', 1, 'box', 'on');

%%% upwelling at PPZ ::
nexttile();
hold('on');
for iC = 1 : 1 : length(IDX_CMEMS)
	plot(xPpz{IDX_CMEMS(iC)}, yPpz{IDX_CMEMS(iC)}, cmems_lineStyles{iC}, ...
	     'color', cmems_color, 'lineWidth', 1, 'displayName', cmems_labels{iC});
end
plot(xPpz{IDX_ECCO}, yPpzcmemsMean,      '-k',                       'lineWidth', 2.5, 'displayName', 'CMEMS Mean');
plot(xPpz{IDX_GREP}, yPpz{IDX_GREP},     '-',  'color', 'r',         'lineWidth', 1.5, 'displayName', 'GREP');
plot(xPpz{IDX_ECCO}, yPpz{IDX_ECCO},     '-',  'color', ecco_color,  'lineWidth', 2,   'displayName', 'ECCO');
plot(xPpz{IDX_ECCO}, yPpzallMean,        '-',  'color', green_color, 'lineWidth', 2.5, 'displayName', 'All-Model Mean');
xline(0, '-k', 'lineWidth', 0.5, 'handleVisibility', 'off');
hold('off');
xlabel('\textbf{Latitude [$^\circ$N]}', 'interpreter', 'latex', 'fontSize', 16);
ylabel('$w$ \textbf{[m d$^{-1}$]}', 'interpreter', 'latex', 'fontSize', 16);
title('\textbf{Kriged Upwelling Velocity at PPZ}', 'interpreter', 'latex', 'fontSize', 18);
legend('location', 'eastOutside', 'interpreter', 'latex', 'fontSize', 12);
set(gca, 'tickLabelInterpreter', 'latex', 'fontSize', 14, 'fontWeight', 'bold', 'lineWidth', 1, 'box', 'on');

set(gcf, 'position', [0, 0, 1400, 900]);
exportgraphics(gcf, [plot_output_basepath 'interpData/upwelling/all_model_transect.pdf'], 'ContentType', 'vector');

%%  figure 2: sensitivity — FOAM and GREP excluded
%   FOAM (GloSea5) is a strong outlier relative to the other four products.
%   GREP is excluded here because it incorporates FOAM in its average.
%   line styles preserve the figure 1 assignments (CGLO '-', GLOR '-.', ORAS ':')
%   so the two figures are visually consistent. ::

%   indices excluding FOAM (index 3 in IDX_CMEMS) and GREP ::
IDX_CMEMS_NF  = [2, 4, 5];    % CGLO, GLOR, ORAS
IDX_INDEP_NF  = [1, 2, 4, 5]; % ECCO, CGLO, GLOR, ORAS

cmems_nf_labels     = {'CGLO', 'GLOR', 'ORAS'};
cmems_nf_lineStyles = {'-', '-.', ':'};   % same assignments as figure 1

%   ensemble means excluding FOAM ::
y100cmemsMeanNF = mean(cat(2, y100{IDX_CMEMS_NF}), 2, 'omitnan');
yPpzcmemsMeanNF = mean(cat(2, yPpz{IDX_CMEMS_NF}), 2, 'omitnan');

y100allMeanNF = mean(cat(2, y100{IDX_INDEP_NF}), 2, 'omitnan');
yPpzallMeanNF = mean(cat(2, yPpz{IDX_INDEP_NF}), 2, 'omitnan');

figure;
tiledlayout(2, 1, 'tileSpacing', 'compact');

%%% upwelling at 100 m (excl. FOAM, GREP) ::
nexttile();
hold('on');
for iC = 1 : 1 : length(IDX_CMEMS_NF)
	plot(x100{IDX_CMEMS_NF(iC)}, y100{IDX_CMEMS_NF(iC)}, cmems_nf_lineStyles{iC}, ...
	     'color', cmems_color, 'lineWidth', 1, 'displayName', cmems_nf_labels{iC});
end
plot(x100{IDX_ECCO}, y100cmemsMeanNF,  '-k',                       'lineWidth', 2.5, 'displayName', 'CMEMS Mean (excl. FOAM)');
plot(x100{IDX_ECCO}, y100{IDX_ECCO},   '-',  'color', ecco_color,  'lineWidth', 2,   'displayName', 'ECCO');
plot(x100{IDX_ECCO}, y100allMeanNF,    '-',  'color', green_color, 'lineWidth', 2.5, 'displayName', 'All-Model Mean (excl. FOAM)');
xline(0, '-k', 'lineWidth', 0.5, 'handleVisibility', 'off');
hold('off');
ylabel('$w$ \textbf{[m d$^{-1}$]}', 'interpreter', 'latex', 'fontSize', 16);
title('\textbf{Kriged Upwelling Velocity at 100 m (excl. FOAM, GREP)}', 'interpreter', 'latex', 'fontSize', 18);
set(gca, 'tickLabelInterpreter', 'latex', 'fontSize', 14, 'fontWeight', 'bold', 'lineWidth', 1, 'box', 'on');

%%% upwelling at PPZ (excl. FOAM, GREP) ::
nexttile();
hold('on');
for iC = 1 : 1 : length(IDX_CMEMS_NF)
	plot(xPpz{IDX_CMEMS_NF(iC)}, yPpz{IDX_CMEMS_NF(iC)}, cmems_nf_lineStyles{iC}, ...
	     'color', cmems_color, 'lineWidth', 1, 'displayName', cmems_nf_labels{iC});
end
plot(xPpz{IDX_ECCO}, yPpzcmemsMeanNF,  '-k',                       'lineWidth', 2.5, 'displayName', 'CMEMS Mean (excl. FOAM)');
plot(xPpz{IDX_ECCO}, yPpz{IDX_ECCO},   '-',  'color', ecco_color,  'lineWidth', 2,   'displayName', 'ECCO');
plot(xPpz{IDX_ECCO}, yPpzallMeanNF,    '-',  'color', green_color, 'lineWidth', 2.5, 'displayName', 'All-Model Mean (excl. FOAM)');
xline(0, '-k', 'lineWidth', 0.5, 'handleVisibility', 'off');
hold('off');
xlabel('\textbf{Latitude [$^\circ$N]}', 'interpreter', 'latex', 'fontSize', 16);
ylabel('$w$ \textbf{[m d$^{-1}$]}', 'interpreter', 'latex', 'fontSize', 16);
title('\textbf{Kriged Upwelling Velocity at PPZ (excl. FOAM, GREP)}', 'interpreter', 'latex', 'fontSize', 18);
legend('location', 'eastOutside', 'interpreter', 'latex', 'fontSize', 12);
set(gca, 'tickLabelInterpreter', 'latex', 'fontSize', 14, 'fontWeight', 'bold', 'lineWidth', 1, 'box', 'on');

set(gcf, 'position', [0, 0, 1400, 900]);
exportgraphics(gcf, [plot_output_basepath 'interpData/upwelling/all_model_transect_noFOAM.pdf'], 'ContentType', 'vector');

%%  end subroutine
