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
    potDensPlot = csaps(potDensDepthProf, potDensProf, P_SPLINE_DENSITY, [], ones(size(potDensProf)));
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
    %[temp10mInit, ~, ~] = objectiveMapping(tempDepthProf, ...
    %    			           tempProf, ...
    %	    			           DELTA_DEPTH, ...
    %	    			           upperTheta, ...
    % 	    			           zeros(size(tempProf)));  % no error from data

    % find temp at delta_depth ::
    [temp10m, ~] = ordinary_kriging([], ...
    				    [], ...
				    tempDepthProf, ...
				    tempProf, ...
    				    [], ...
	    			    [], ...
	    			    [DELTA_DEPTH], ...
    				   'FitHp', true, ...
    				   'Variance', 1.0, ...
	    			   'NoiseVariance', 0.1, ...
    				   'GpVerticalScaleM', 100, ...
    				   'RandomState', 7);  % S3: seed for reproducible hyperparameter optimisation

    %   print out new ::
    % disp(['Station #', num2str(sn), '- Objective mapping: T(10m) = ', num2str(round(temp10mInit, 3)), ' | Ordinary Kriging: T(10m) = ', num2str(round(temp10m, 3))]);

    %   find threshold value ::
    delta_t = temp10m - DELTA_T;                               
    
    %   best fit spline ::
    tempProfFit = csaps(tempDepthProf, tempProf - delta_t, P_SPLINE_DENSITY, [], ones(size(tempProf)));
    
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
    potDensFit = csaps(potDensDepthProf, potDensProf - mldPPotDens, P_SPLINE_DENSITY, [], ones(size(potDensProf)));
    
    %   find mld depth given potDens std :: 
    pZeros = fnzeros(potDensFit); % change to splineroots
    mld_JAK(iStat) = min(pZeros(pZeros > MLD_AVE_DEPTH), [], 'all'); 
    
    % ---------- BW method ------------  
    %   get surface value ::
    potDensSurf = mean(potDensProf(potDensDepthProf <= MLD_AVE_DEPTH), 'all');

    %   calculate sigma-theta MLD (B&W) ::
    for iDelta = 1 : 1 : length(DELTA_POTDENS)
        mldBWPotDens = potDensSurf + DELTA_POTDENS(iDelta);
        potDensFit = csaps(potDensDepthProf, potDensProf - mldBWPotDens, P_SPLINE_DENSITY, [], ones(size(potDensProf)));
        BWZeros = fnzeros(potDensFit); 
        if isempty(BWZeros(BWZeros > MLD_AVE_DEPTH))
            mld_BW(iStat, iDelta) = min(BWZeros, [], 'all'); 
        else 
            mld_BW(iStat, iDelta) = min(BWZeros(BWZeros > MLD_AVE_DEPTH), [], 'all'); 
        end
    end
    
    % see plotMld.m for diagnostic figures

end

%% end subroutine
