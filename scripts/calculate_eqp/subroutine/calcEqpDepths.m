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
   
    %   valid data for ratio fit (remove NaN; guard th234 > 0 and u238 > 0) ::
    validIdx   = ~isnan(dataSn.th234) & ~isnan(dataSn.u238) & ...
                 ~isnan(dataSn.uncertTh234) & ~isnan(dataSn.uncertU238) & ...
                 (dataSn.th234 > 0) & (dataSn.u238 > 0);
    depthV     = dataSn.depth(validIdx);
    th234V     = dataSn.th234(validIdx);
    u238V      = dataSn.u238(validIdx);

    %   ratio, uncertainty (delta method for division), and inverse-variance weights ::
    ratioV     = th234V ./ u238V;
    diff       = ratioV - 1;
    sigmaRatio = ratioV .* sqrt((dataSn.uncertTh234(validIdx) ./ th234V).^2 + ...
                                (dataSn.uncertU238(validIdx)  ./ u238V ).^2);
    wts_eqp    = 1 ./ (sigmaRatio .^ 2);

    %   csaps spline fits with inverse-variance weights ::
    % fitDiff    = fit(dataSn.depth, diff,       'smoothingspline');  % old: different toolbox, no weights, version-dependent p
    % fitDiffAbs = fit(dataSn.depth, abs(diff),  'smoothingspline');
    fitDiff    = csaps(depthV, diff,       P_SPLINE_GRADIENT, [], wts_eqp);
    fitDiffAbs = csaps(depthV, abs(diff),  P_SPLINE_GRADIENT, [], wts_eqp);

    %   find roots ::
    % roots{iStat} = fnzeros(fitDiff.p);  % old: .p extracted pp struct from cfit object
    roots{iStat} = fnzeros(fitDiff);
    rootStat     = roots{iStat}(1, :);

    %   find minimum of |ratio - 1| ::
    % mins{iStat} = fminbnd(fitDiffAbs, X0, ez + TOL);  % old: cfit object callable directly
    mins{iStat} = fminbnd(@(z) ppval(fitDiffAbs, z), X0, ez + TOL);
    minStat     = mins{iStat}(1, :);

    %   find final points ::
    W_EZ = 0.5;
    W_MLD = 0.5;
    possVals = [rootStat, minStat]';
    distVals = (W_EZ .* abs(ez - possVals)) + (W_MLD .* abs(mldVal - possVals));
    [~, closestIndex] = min(distVals);
    rootsFt(iStat) = possVals(closestIndex);
    
end

%   save data ::
%%% excel ::
gp15_eqp = array2table([gp15_stations.stationNo, rootsFt]);
gp15_eqp.Properties.VariableNames = {'stationNo', 'eqp'};
writetable(gp15_eqp, [sim_output_basepath 'calcEqp/gp15_eqp.xlsx'])

%%% mat ::
save([sim_output_basepath 'calcEqp/gp15_eqp.mat'], 'gp15_eqp');

%%  end subroutine
