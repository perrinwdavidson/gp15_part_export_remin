%%  calculate final mld
%   preallocate ::
mld_JAK = NaN(NUMSTAT, 1);
mld_dBM = NaN(NUMSTAT, 1);
mld_BW = NaN(NUMSTAT, length(DELTA_POTDENS));

%   set correlation function and values ::
% corFun = 'sexp';
% lowerTheta = 0;
upperTheta = 10;

%   run through stations ::
for iStat = 1 : 1 : NUMSTAT  
    
    %   find station ::
    sn = gp15_stations.stationNo(iStat); 
    
    %   find cast ::
    statCast = gp15_stations.castNo(iStat); 
    
    %   get mld guess ::
    mldGuess = mldCalcs.mld_JAK(iStat); 

    %   only use first downcast ::
    fullDepthPlot = gp15_ctd.ctd_depth(gp15_ctd.stationNo == sn & ...
                                       gp15_ctd.castNo == statCast);
    iDepthPlot = find(fullDepthPlot >= mldGuess + ADD_DEPTH);
    if isempty(iDepthPlot)
        endDepthPlot = length(fullDepthPlot);
    else
        endDepthPlot = iDepthPlot(1);
    end

    %   make depth appropriate array ::
    iCtdStatCast = find(gp15_ctd.stationNo == sn & gp15_ctd.castNo == statCast);
    ctdStatCast = gp15_ctd(iCtdStatCast, :);
    ctdStatCastDepth = ctdStatCast(1 : endDepthPlot, :);

    %   find profiles ::
    potDensProf = ctdStatCastDepth.ctd_potentialDensity;
    potDensDepthProf = ctdStatCastDepth.ctd_depth;
    
    %   fit for plotting ::
    potDensPlot = csaps(potDensDepthProf, potDensProf); 
    potDensPlotDepth = 0 : 0.1 : max(potDensDepthProf);
    potDensPlot = ppval(potDensPlot, potDensPlotDepth); 
    
    % ---------- dBM method ------------
    %   find profiles ::
    tempProf = gp15_ctd.ctd_temperature(gp15_ctd.stationNo == sn & ...
                                        gp15_ctd.castNo == statCast);
    tempDepthProf = gp15_ctd.ctd_depth(gp15_ctd.stationNo == sn & ...
                                       gp15_ctd.castNo == statCast);

    %   find temp at delta_depth ::
    %   n.b.: we assume that all points have error (including the value given), unlike kriging below.
    % [temp10m, ~, ~, ~, ~] = krigeData(tempDepthProf, ...
    %                                   tempProf, ...
    %                                   DELTA_DEPTH, ...
    %                                   corFun, ...
    %                                   lowerTheta, ...
    %                                   upperTheta);   

    %   find temp at delta_depth ::
    [temp10m, ~, ~] = objectiveMapping(tempDepthProf, ...
	    			       tempProf, ...
	    			       DELTA_DEPTH, ...
	    			       upperTheta, ...
	    			       zeros(size(tempProf)));  % no error from data

    %   find threshold value ::
    delta_t = temp10m - DELTA_T;                               
    
    %   best fit spline ::
    tempProfFit = csaps(tempDepthProf, tempProf - delta_t); 
    
    %   find depth ::
    dBMZeros = fnzeros(tempProfFit);
    mld_dBM(iStat) = min(dBMZeros(dBMZeros > MLD_AVE_DEPTH), [], 'all'); 
    
    % ---------- Pickart method ------------
    %   subjective mld ::
    mldSubjective = mldCalcs.mld_JAK(iStat);
    
    %   mld arrays of values below this depth ::
    potDensMld = potDensProf(potDensDepthProf <= mldSubjective);
    
    %   find statistics of mld values ::
    statsMld.meanMld(iStat) = mean(potDensMld, 'all');
    statsMld.stdMld(iStat) = std(potDensMld, 0, 'all');                     % 0 means 1/N-1 (sample variance used)
    
    %   find mld depth envelope given std ::
    mldPPotDens = statsMld.meanMld(iStat) + (NUM_STD * statsMld.stdMld(iStat));
    
    %   find spline best fit ::
    potDensFit = csaps(potDensDepthProf, potDensProf - mldPPotDens); 
    
    %   find mld depth given potDens std :: 
    pZeros = fnzeros(potDensFit); % change to splineroots
    mld_JAK(iStat) = min(pZeros(pZeros > MLD_AVE_DEPTH), [], 'all'); 
    
    % ---------- BW method ------------  
    %   get surface value ::
    potDensSurf = mean(potDensProf(potDensDepthProf <= MLD_AVE_DEPTH), 'all');

    %   calculate sigma-theta MLD (B&W) ::
    for iDelta = 1 : 1 : length(DELTA_POTDENS)
        mldBWPotDens = potDensSurf + DELTA_POTDENS(iDelta);
        potDensFit = csaps(potDensDepthProf, potDensProf - mldBWPotDens); 
        BWZeros = fnzeros(potDensFit); 
        if isempty(BWZeros(BWZeros > MLD_AVE_DEPTH))
            mld_BW(iStat, iDelta) = min(BWZeros, [], 'all'); 
        else 
            mld_BW(iStat, iDelta) = min(BWZeros(BWZeros > MLD_AVE_DEPTH), [], 'all'); 
        end
    end
    
    %   plot profile ::
    if strcmp(plotting, 'yes')
        
        figure; 

        hold on

        scatter(potDensProf, potDensDepthProf, 's', ...
                'lineWidth', 1, 'markerEdgeColor', 'k', ...
                'markerFaceColor', 'white');
        plot(potDensPlot, potDensPlotDepth, '--k', ...
             'markerEdgeColor', 'k', 'markerFaceColor', 'w', ...
             'lineWidth', 1, 'markerSize', 5);

        yline(mld_JAK(iStat), '-k', 'lineWidth', 2)

        yline(mld_dBM(iStat), '-k', 'lineWidth', 1)
        yline(mld_BW(iStat, 1), '--k', 'LineWidth', 1)
        yline(mld_BW(iStat, 2), '-.k', 'LineWidth', 1)
        yline(mld_BW(iStat, 3), ':k', 'LineWidth', 1)
        xline(statsMld.meanMld(iStat) + (NUM_STD * statsMld.stdMld(iStat)), ...
              ':k', 'lineWidth', 0.5)
        xline(statsMld.meanMld(iStat) - (NUM_STD * statsMld.stdMld(iStat)), ...
              ':k', 'lineWidth', 0.5)

        hold off

        box on

        xlabel('\textbf{Potential Density  (kg m$^{-3}$)}', ...
               'interpreter', 'latex', 'fontSize', 20); 
        ylabel('\textbf{Depth (m)}', 'interpreter', 'latex', ...
               'fontSize', 20); 

        ylim([0  80])

        set(gca, 'yDir', 'reverse', 'tickLabelInterpreter', ...
            'latex', 'fontSize', 16, 'fontWeight', 'bold', ...
            'lineWidth', 1);

        set(gcf, 'units', 'inches', 'position', [2, 2, 5, 10], ...
            'paperUnits', 'inches', 'paperSize', [5, 10]);

        exportgraphics(gcf, ...
                       [plot_output_basepath 'calcMld/mldPlots/mldcalc' num2str(sn) '.png'], ...
                       'resolution', 300); 
                   
    end

end

%% end subroutine
