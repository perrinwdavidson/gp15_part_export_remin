%%  plot equatorial stations
%   set whether or not you want to plot ::
plotting = 'yes';

%   set max depth ::
maxDepth = 10000; 

%   make output array ::
gp15_grad = gp15_obs(:, {'stationNo', 'longitude', 'latitude', 'depth'}); 
gp15_grad.vertGrad = NaN(size(gp15_obs, 1), 1);
gp15_grad.vertGradError = NaN(size(gp15_obs, 1), 1);

%   loop through all stations ::
for iStat = 1 : 1 : NUMSTAT
    
    % find station:
    sn = gp15_stations.stationNo(iStat);
    
    %   get station number ::
    statNo = gp15_stations.stationNo(iStat); 
    
    % find ppz:
    ez = gp15_ppz.ppzDepth(gp15_ppz.stationNo == sn);

    % find mld ::
    mldVal = gp15_mld.("MLD JAK")(iStat);
    
    % roots ::
    rootsFt = gp15_eqp.eqp(gp15_eqp.stationNo == sn);
        
    %   get data ::
    x = gp15_obs.depth(gp15_obs.stationNo == statNo &  gp15_obs.depth <= maxDepth); 
    yth = gp15_obs.th234(gp15_obs.stationNo == statNo &  gp15_obs.depth <= maxDepth); 
    nth = gp15_obs.uncertTh234(gp15_obs.stationNo == statNo &  gp15_obs.depth <= maxDepth); 
    
    %   calculate number of depths ::
    n = length(x); 
    
    %   calculate gradient error ::
    for iDepth = 1 : 1 : n
        
        %   calculate for forward difference ::
        if iDepth == 1
            
            %   get thorim values ::
            th1 = yth(iDepth); 
            th2 = yth(iDepth + 1); 
            
            %   get thorium error values ::
            deltaTh1 = nth(iDepth); 
            deltaTh2 = nth(iDepth + 1); 
            
            %   get depth values ::
            z1 = x(iDepth); 
            z2 = x(iDepth + 1);
            
            %   calculate B ::
            B = abs(z2 - z1); 
            
            %   calculate gradient ::
            gradTh(iDepth) = (th2 - th1) / (z2 - z1); 
            
            %   calculate error ::
            deltaGradTh(iDepth) = sqrt(((deltaTh1 / B) ^ 2) ...
                                       + ((deltaTh2 / B) ^ 2)); 
            
        %   calculate for backward difference ::
        elseif iDepth == n
            
            %   get thorim values ::
            th1 = yth(iDepth - 1); 
            th2 = yth(iDepth); 
            
            %   get thorium error ::
            deltaTh1 = nth(iDepth - 1); 
            deltaTh2 = nth(iDepth); 
            
            %   get depth values ::
            z1 = x(iDepth - 1); 
            z2 = x(iDepth);
            
            %   calculate B ::
            B = abs(z2 - z1); 
            
            %   calculate gradient ::
            gradTh(iDepth) = (th2 - th1) / (z2 - z1); 
            
            %   calculate error ::
            deltaGradTh(iDepth) = sqrt(((deltaTh1 / B) ^ 2) ...
                                       + ((deltaTh2 / B) ^ 2)); 
            
        %   calculate for centered difference ::
        else
            
            %   get thorim values ::
            th1 = yth(iDepth - 1); 
            th2 = yth(iDepth + 1); 
            
            %   get thorium error ::
            deltaTh1 = nth(iDepth - 1); 
            deltaTh2 = nth(iDepth + 1); 
            
            %   get depth values ::
            z1 = x(iDepth - 1); 
            z2 = x(iDepth + 1);
            
            %   calculate B ::
            B = abs(z2 - z1); 
            
            %   calculate gradient ::
            gradTh(iDepth) = (th2 - th1) / (z2 - z1); 
            
            %   calculate error ::
            deltaGradTh(iDepth) = sqrt(((deltaTh1 / B) ^ 2) ...
                                       + ((deltaTh2 / B) ^ 2)); 
            
        end
        
    end
    deltaGradTh = deltaGradTh .^ 2;  % make variance 
    
    %   smoothing spline ::
    % weights = 1 ./ deltaGradTh;
    % gradThHat = csaps(x, gradTh, [], x, weights);
    
    %    other error versions ::
    % deltaGradThHat = abs(gradThHat - gradTh');  % <- best option
    % deltaGradThHat = sqrt((deltaGradThHat .^ 2) + (deltaGradTh .^ 2));

    %   moving mean ::
    lz = 25; 
    n = length(gradTh); 
    gradThHat = NaN(size(gradTh)); 
    deltaGradThHat = NaN(size(gradTh)); 
    for i = 1 : 1 : n
	ix = x(i); 
	xStart = ix - (lz / 2); 
	xEnd = ix + (lz / 2); 
	gradThGet = gradTh((x >= xStart) & (x <= xEnd));
	deltaGradThGet = deltaGradTh((x >= xStart) & (x <= xEnd));
	gradThHat(i) = mean(gradThGet); 
	deltaGradThHat(i) = mean(deltaGradThGet); 
    end
    deltaGradThHat = sqrt(deltaGradThHat);  % make standard deviation

    %   plotting arrays ::
    if strcmp(plotting, 'yes')

        figure; 

        hold on

        plot(gradTh, x, '-ok', ...
             'markerEdgeColor', 'k', ...
             'markerFaceColor', 'k', ...
             'lineWidth', 0.25, ...
             'markerSize', 5);
        plot(gradThHat, x, '--', 'color', 'k', 'lineWidth', 1);
        scatter(deltaGradThHat, x, 'o', ...
               'lineWidth', 1, ...
               'markerEdgeColor', 'k', ...
               'markerFaceColor', 'white');

        hold off

        box on

        xlabel('\textbf{Gradient (dpm L$^{-1}$ m$^{-1}$)}', ...
               'interpreter', 'latex', 'fontSize', 20); 
        ylabel('\textbf{Depth (m)}', 'interpreter', 'latex', ...
               'fontSize', 20); 

        ylim([0  200]);
        xlim([min(gradTh, [], 'all'), max(gradTh, [], 'all')]);

        set(gca, 'yDir', 'reverse', 'tickLabelInterpreter', ...
            'latex', 'fontSize', 16, 'fontWeight', 'bold', ...
            'lineWidth', 1);

        set(gcf, 'units', 'inches', 'position', [2, 2, 5, 10], ...
            'paperUnits', 'inches', 'paperSize', [5, 10]);

        exportgraphics(gcf, ...
                       [plot_output_basepath 'calcVertGrad/gradPlots/grad' num2str(statNo) '.png'], ...
                       'resolution', 300); 

    end

    %   remove values outside of gradient:
    gradThHat(x < mldVal | x > rootsFt | x > ez) = 0;
    deltaGradThHat(x < mldVal | x > rootsFt | x > ez) = 0;
    
    %   make final arrays
    %%% get indices of where data is ::
    statIdx = find(gp15_obs.stationNo == statNo & gp15_obs.depth <= maxDepth); 

    %%% put data in ::
    gp15_grad.vertGrad(statIdx) = gradThHat; 
    gp15_grad.vertGradError(statIdx) = deltaGradThHat; 
            
    %   clear arrays ::
    clear('gradTh', 'deltaGradTh'); 
            
end
close('all'); 

%%  end subroutine
