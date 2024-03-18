%   do you want to plot? yes or no. 
iplot = 'yes'; 

%   initialize ::
%%% make pn lsf ::
len_data = size(gp15_obs, 1); 
gp15_obs.pnLarge = NaN(len_data, 1);
gp15_obs.uncertPnLarge = NaN(len_data, 1);
gp15_obs.pnTh234RatioLarge = NaN(len_data, 1);
gp15_obs.uncertPnTh234RatioLarge = NaN(len_data, 1);

%%% make pn ssf ::
gp15_obs.pnSmall = NaN(len_data, 1);
gp15_obs.uncertPnSmall = NaN(len_data, 1);
gp15_obs.pnTh234RatioSmall = NaN(len_data, 1);
gp15_obs.uncertPnTh234RatioSmall = NaN(len_data, 1);

%   loop through all stations ::
for iStat = 1 : 1 : NUMSTAT

	%   get station number ::
	statNo = gp15_stations.stationNo(iStat);

	%   get station data ::
	statData = gp15_obs(gp15_obs.stationNo == statNo, :);
	statDataRadioPn = gp15_part(gp15_part.stationNo == statNo, :);
	statDataPn = gp15_particles(gp15_particles.station_num == statNo, :);  

	%    get depths ::
	xq = statData.depth;

	%    get sample data ::
	XPnLarge = table2array(sortrows(rmmissing(statDataPn(:, {'depth_m', 'PN_LPT_uM', 'PN_LPT_uM_std'})), 'depth_m'));
	XPnSmall = table2array(sortrows(rmmissing(statDataPn(:, {'depth_m', 'PN_SPT_uM', 'PN_SPT_uM_std'})), 'depth_m'));

	%    deal with empty arrays ::
	if isempty(XPnLarge)
		XPnLarge = NaN(size(statData, 1), 3);
	end
	if isempty(XPnSmall)
		XPnSmall = NaN(size(statData, 1), 3);
	end

	%    interpolate ::
	[pnLarge, pnLargeStd] = interp1dError(XPnLarge(:, 1), XPnLarge(:, 2), XPnLarge(:, 3), xq); 
	[pnSmall, pnSmallStd] = interp1dError(XPnSmall(:, 1), XPnSmall(:, 2), XPnSmall(:, 3), xq); 
	
	%    make ratio data :: 
	pnTh234RatioLarge = pnLarge ./ statData.th234PocLarge;
	pnTh234RatioLargeStd = abs(pnTh234RatioLarge) .* sqrt(((pnLargeStd ./ pnLarge) .^ 2) + ((statData.uncertTh234PocLarge ./ statData.th234PocLarge) .^ 2));
	pnTh234RatioSmall = pnSmall ./ statData.th234PocSmall;
	pnTh234RatioSmallStd = abs(pnTh234RatioSmall) .* sqrt(((pnSmallStd ./ pnSmall) .^ 2) + ((statData.uncertTh234PocSmall ./ statData.th234PocSmall) .^ 2));

        %    plot ::
        if strcmp(iplot, 'yes') 
	
            figure;
	    tl = tiledlayout(2, 4, 'tileSpacing', 'compact');
       
	    % u238 ::
	    nexttile();

            hold('on');
            plot(statData.u238, statData.depth, '-ok', ...
                 'markerEdgeColor', 'k', 'markerFaceColor', 'k', ...
                 'lineWidth', 1, 'markerSize', 5);
            scatter(statData.uncertU238, statData.depth, 's', ...
                    'lineWidth', 1, 'markerEdgeColor', 'k', ...
                    'markerFaceColor', 'white');
            hold('off');
            
            box('on');
            
            xlabel('\textbf{$^{238}$U (dpm L$^{-1}$)}', ...
                   'interpreter', 'latex', 'fontSize', 20); 
            ylabel('\textbf{Depth (m)}', 'interpreter', 'latex', ...
                   'fontSize', 20); 
               
            ylim([0 1000]);
              
            set(gca, 'yDir', 'reverse', 'tickLabelInterpreter', ...
                'latex', 'fontSize', 16, 'fontWeight', 'bold', ...
                'lineWidth', 1);

	    % th234 ::
	    nexttile();

            hold('on');
            plot(statData.th234, statData.depth, '-ok', ...
                 'markerEdgeColor', 'k', 'markerFaceColor', 'k', ...
                 'lineWidth', 1, 'markerSize', 5);
            scatter(statData.uncertTh234, statData.depth, 's', ...
                    'lineWidth', 1, 'markerEdgeColor', 'k', ...
                    'markerFaceColor', 'white');
            hold('off');
            
            box('on');
            
            xlabel('\textbf{Total $^{234}$Th (dpm L$^{-1}$)}', ...
                   'interpreter', 'latex', 'fontSize', 20); 
            
	    ylim([0 1000]);
              
            set(gca, 'yDir', 'reverse', 'tickLabelInterpreter', ...
                'latex', 'fontSize', 16, 'fontWeight', 'bold', ...
                'lineWidth', 1);

	    % th234 ssf ::
	    nexttile();

            hold('on');
            plot(statData.th234PocSmall, statData.depth, '-ok', ...
                 'markerEdgeColor', 'k', 'markerFaceColor', 'k', ...
                 'lineWidth', 1, 'markerSize', 5);
            scatter(statData.uncertTh234PocSmall, statData.depth, 's', ...
                    'lineWidth', 1, 'markerEdgeColor', 'k', ...
                    'markerFaceColor', 'white');
            hold('off');
            
            box('on');
            
            xlabel('\textbf{SSF Particulate $^{234}$Th (dpm L$^{-1}$)}', ...
                   'interpreter', 'latex', 'fontSize', 20); 
               
            ylim([0 1000]);
              
            set(gca, 'yDir', 'reverse', 'tickLabelInterpreter', ...
                'latex', 'fontSize', 16, 'fontWeight', 'bold', ...
                'lineWidth', 1);

	    % th234 lsf ::
	    nexttile();

            hold('on');
            plot(statData.th234PocLarge, statData.depth, '-ok', ...
                 'markerEdgeColor', 'k', 'markerFaceColor', 'k', ...
                 'lineWidth', 1, 'markerSize', 5);
            scatter(statData.uncertTh234PocLarge, statData.depth, 's', ...
                    'lineWidth', 1, 'markerEdgeColor', 'k', ...
                    'markerFaceColor', 'white');
            hold('off');
            
            box('on');
            
            xlabel('\textbf{LSF Particulate $^{234}$Th (dpm L$^{-1}$)}', ...
                   'interpreter', 'latex', 'fontSize', 20); 
               
            ylim([0 1000]);
              
            set(gca, 'yDir', 'reverse', 'tickLabelInterpreter', ...
                'latex', 'fontSize', 16, 'fontWeight', 'bold', ...
                'lineWidth', 1);

	    % ssf pn ::
	    nexttile();

            hold('on');
            plot(pnSmall, statData.depth, '-ok', ...
                 'markerEdgeColor', 'k', 'markerFaceColor', 'k', ...
                 'lineWidth', 1, 'markerSize', 5);
            scatter(pnSmallStd, statData.depth, 's', ...
                    'lineWidth', 1, 'markerEdgeColor', 'k', ...
                    'markerFaceColor', 'white');
            hold('off');
            
            box('on');
            
            xlabel('\textbf{SSF PN ($\mu$M)}', ...
                   'interpreter', 'latex', 'fontSize', 20); 
            ylabel('\textbf{Depth (m)}', 'interpreter', 'latex', ...
                   'fontSize', 20); 
               
            ylim([0 1000]);
              
            set(gca, 'yDir', 'reverse', 'tickLabelInterpreter', ...
                'latex', 'fontSize', 16, 'fontWeight', 'bold', ...
                'lineWidth', 1);

	    % lsf pn ::
	    nexttile();

            hold('on');
            plot(pnLarge, statData.depth, '-ok', ...
                 'markerEdgeColor', 'k', 'markerFaceColor', 'k', ...
                 'lineWidth', 1, 'markerSize', 5);
            scatter(pnLargeStd, statData.depth, 's', ...
                    'lineWidth', 1, 'markerEdgeColor', 'k', ...
                    'markerFaceColor', 'white');
            hold('off');
            
            box('on');
            
            xlabel('\textbf{LSF PN ($\mu$M)}', ...
                   'interpreter', 'latex', 'fontSize', 20); 
               
            ylim([0 1000]);
              
            set(gca, 'yDir', 'reverse', 'tickLabelInterpreter', ...
                'latex', 'fontSize', 16, 'fontWeight', 'bold', ...
                'lineWidth', 1);

	    % ssf pn:234th ::
	    nexttile();

            hold('on');
            plot(pnTh234RatioSmall, statData.depth, '-ok', ...
                 'markerEdgeColor', 'k', 'markerFaceColor', 'k', ...
                 'lineWidth', 1, 'markerSize', 5);
            scatter(pnTh234RatioSmallStd, statData.depth, 's', ...
                    'lineWidth', 1, 'markerEdgeColor', 'k', ...
                    'markerFaceColor', 'white');
            hold('off');
            
            box('on');
            
            xlabel('\textbf{SSF PN:$^{234}$Th ($\mu$mol dpm$^{-1}$)}', ...
                   'interpreter', 'latex', 'fontSize', 20); 
               
            ylim([0 1000]);
              
            set(gca, 'yDir', 'reverse', 'tickLabelInterpreter', ...
                'latex', 'fontSize', 16, 'fontWeight', 'bold', ...
                'lineWidth', 1);

	    % lsf pn:234th ::
	    nexttile();

            hold('on');
            plot(pnTh234RatioLarge, statData.depth, '-ok', ...
                 'markerEdgeColor', 'k', 'markerFaceColor', 'k', ...
                 'lineWidth', 1, 'markerSize', 5);
            scatter(pnTh234RatioLargeStd, statData.depth, 's', ...
                    'lineWidth', 1, 'markerEdgeColor', 'k', ...
                    'markerFaceColor', 'white');
            hold('off');
            
            box('on');
            
            xlabel('\textbf{LSF PN:$^{234}$Th ($\mu$mol dpm$^{-1}$)}', ...
                   'interpreter', 'latex', 'fontSize', 20); 
               
            ylim([0 1000]);
              
            set(gca, 'yDir', 'reverse', 'tickLabelInterpreter', ...
                'latex', 'fontSize', 16, 'fontWeight', 'bold', ...
                'lineWidth', 1);

            title(tl, ['\textbf{Station ' num2str(statNo) '}'], ...
                  'interpreter', 'latex', 'fontSize', 24); 

            set(gcf, 'units', 'inches', 'position', [0, 0, 20, 15], ...
                'paperUnits', 'inches', 'paperSize', [20, 15]);
           
            exportgraphics(gcf, [plot_output_basepath, 'doQc/stationPlots/pn/station', num2str(statNo), '.png'], ...
                           'resolution', 300); 
            
        end
        
	%   put into data array ::
	statData.pnLarge = pnLarge;
	statData.uncertPnLarge = pnLargeStd;
	statData.pnTh234RatioLarge = pnTh234RatioLarge;
	statData.uncertPnTh234RatioLarge = pnTh234RatioLargeStd;
	statData.pnSmall = pnSmall;
	statData.uncertPnSmall = pnSmallStd;
	statData.pnTh234RatioSmall = pnTh234RatioSmall;
	statData.uncertPnTh234RatioSmall = pnTh234RatioSmallStd;

	%   put into storage array ::
	if iStat == 1
		gp15_obs_final = statData;
	else
		gp15_obs_final = [gp15_obs_final; statData];
	end

end

%%  hand off data
gp15_obs = gp15_obs_final;

%%  clean up
clear('gp15_obs_final');
close('all'); 

%%  end subsubroutine
