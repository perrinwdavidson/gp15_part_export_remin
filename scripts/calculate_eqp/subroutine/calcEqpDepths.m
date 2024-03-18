%%  calculate and plot equilibrium points (eqp)
%   preallocate ::
roots = cell(NUMSTAT, 1);
mins = cell(NUMSTAT, 1);
rootsFt = NaN(NUMSTAT, 1);

for iStat = 1 : 1 : NUMSTAT  

    %   find station ::
    sn = gp15_stations.stationNo(iStat);

    %   find ppz ::
    ez = gp15_ppz.ppzDepth(gp15_ppz.stationNo == sn);

    %   find mld ::
    mldVal = gp15_mld.("MLD JAK")(iStat);    

    %   find profile ::
    dataSn = gp15_obs(gp15_obs.stationNo == sn & gp15_obs.depth <= BTM_DEPTH, :);
   
    %   difference ::    
    diff = (dataSn.th234 ./ dataSn.u238) - 1;
    fitDiff = fit(dataSn.depth, diff, 'smoothingspline');
    fitDiffAbs = fit(dataSn.depth, abs(diff), 'smoothingspline'); 
    
    %   find roots ::
    roots{iStat} = fnzeros(fitDiff.p);
    rootStat = roots{iStat}(1, :);

    %   find minimum points ::
    mins{iStat} = fminbnd(fitDiffAbs, X0, ez + TOL);
    minStat = mins{iStat}(1, :);

    %   find final points ::
    W_EZ = 0.5;
    W_MLD = 0.5;
    possVals = [rootStat, minStat]';
    distVals = (W_EZ .* abs(ez - possVals)) + (W_MLD .* abs(mldVal - possVals));
    [~, closestIndex] = min(distVals);
    rootsFt(iStat) = possVals(closestIndex);
    
    %   plotting arrays ::
    if strcmp(plotting, 'yes')
        
        figure; 

        hold on

        plot(dataSn.th234, dataSn.depth, '-ok', ...
                'lineWidth', 1, 'markerEdgeColor', 'k', ...
                'markerFaceColor', 'white');
        plot(dataSn.u238, dataSn.depth, '-sk', ...
             'markerEdgeColor', 'k', 'markerFaceColor', 'w', ...
             'lineWidth', 1, 'markerSize', 5);

        yline(rootsFt(iStat), '-k', 'lineWidth', 2)
        yline(mldVal, '--k', 'lineWidth', 1)
        yline(ez, '-.k', 'LineWidth', 1)

        hold off

        box on

        xlabel('\textbf{Activity  (dpm L$^{-1}$)}', ...
               'interpreter', 'latex', 'fontSize', 20); 
        ylabel('\textbf{Depth (m)}', 'interpreter', 'latex', ...
               'fontSize', 20); 

        ylim([0  200])
        xlim([1 3])

        set(gca, 'yDir', 'reverse', 'tickLabelInterpreter', ...
            'latex', 'fontSize', 16, 'fontWeight', 'bold', ...
            'lineWidth', 1);

        set(gcf, 'units', 'inches', 'position', [2, 2, 5, 10], ...
            'paperUnits', 'inches', 'paperSize', [5, 10]);

        exportgraphics(gcf, ...
                       [plot_output_basepath 'calcEqp/eqpPlots/eqpcalc' num2str(sn) '.png'], ...
                       'resolution', 300); 
                   
    end

end

%   save data ::
%%% excel ::
gp15_eqp = array2table([gp15_stations.stationNo, rootsFt]);
gp15_eqp.Properties.VariableNames = {'stationNo', 'eqp'};
writetable(gp15_eqp, [sim_output_basepath 'calcEqp/gp15_eqp.xlsx'])

%%% mat ::
save([sim_output_basepath 'calcEqp/gp15_eqp.mat'], 'gp15_eqp');

%%  clean up
close('all');
%%  end subroutine
