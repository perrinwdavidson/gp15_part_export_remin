%%  temporal average (NOT SAVED -- DONE AFTER INTERPOLATION)
%   get indices ::
idxLon = find(min(abs(u_mercator.longitude - mean(gp15_stations.longitude + 360))) == abs(u_mercator.longitude - mean(gp15_stations.longitude + 360)));
idxDepth = find(min(abs(u_mercator.depth - 100)) == abs(u_mercator.depth - 100));
idxTime = find(min(abs(u_mercator.time - mean(gp15_stations.date))) == abs(u_mercator.time - mean(gp15_stations.date)));

%   get data ::
wAveDay = squeeze(mean(w(idxLon, :, 1:idxDepth, idxTime), 3))' * SEC2DAY; 

%   smooth again ::
wAveDay = movmean(wAveDay, degreeAve / spaceResolution, 1); 

%%  plot 
%   find latitude ::
wLat = u_mercator.latitude;

%   plot ::
figure;
plot(wLat, wAveDay, '-x', 'lineWidth', 1, 'color', 'k');
title(['\textbf{Spatially and Temporally Averaged PMT Vertical Advection Velocity at 100 m (' upper(dataProduct) ')}'], 'fontSize', 18, 'fontWeight', 'bold', 'interpreter', 'latex');
xlabel('\textbf{Latitude (Degrees North)}', 'fontSize', 18, 'interpreter', 'latex');
ylabel('\textbf{Vertical Velocity (m day$^{-1}$)}', 'fontSize', 18, 'interpreter', 'latex');
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 20, 8], 'PaperUnits', 'Inches', 'PaperSize', [20, 8]);
exportgraphics(gcf, [plot_output_basepath 'calcUpwell/w_' dataProduct '.png'], 'resolution', 300);  

%%  end subroutine
