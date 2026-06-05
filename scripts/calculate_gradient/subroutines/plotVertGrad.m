%% plotVertGrad - diagnostic plots of vertical Th-234 gradient profiles
%  standalone; run after calcVertGrad.m has saved gp15_grad.mat.
%  recomputes raw finite differences from gp15_obs for comparison against
%  the saved smoothed + masked gradient.
%--------------------------------------------------------------------------
setupModel;
setModelCoefficients;

%%  load saved outputs
load([pro_output_basepath 'interpData/observations/gp15_obs.mat'], 'gp15_obs');
load([sim_output_basepath 'calcVertGrad/gp15_grad.mat'], 'gp15_grad');
load([pro_output_basepath 'doQc/ppz/gp15_ppz.mat'], 'gp15_ppz');
load([sim_output_basepath 'calcMld/gp15_mld.mat'], 'gp15_mld');
load([sim_output_basepath 'calcEqp/gp15_eqp.mat'], 'gp15_eqp');
load([pro_output_basepath 'doQc/stations/gp15_stations.mat'], 'gp15_stations', 'NUMSTAT');

maxDepth = 10000;

%%  loop over stations
for iStat = 1 : 1 : NUMSTAT

    sn = gp15_stations.stationNo(iStat);

    x   = gp15_obs.depth(gp15_obs.stationNo == sn & gp15_obs.depth <= maxDepth);
    yth = gp15_obs.th234(gp15_obs.stationNo == sn & gp15_obs.depth <= maxDepth);
    nth = gp15_obs.uncertTh234(gp15_obs.stationNo == sn & gp15_obs.depth <= maxDepth);
    n   = length(x);

    %   raw finite differences (deterministic from gp15_obs; identical to doVertGradCalc) ::
    gradTh      = zeros(n, 1);
    deltaGradTh = zeros(n, 1);
    for iDepth = 1 : 1 : n
        if iDepth == 1
            th1 = yth(1);   th2 = yth(2);
            dt1 = nth(1);   dt2 = nth(2);
            z1  = x(1);     z2  = x(2);
        elseif iDepth == n
            th1 = yth(n-1); th2 = yth(n);
            dt1 = nth(n-1); dt2 = nth(n);
            z1  = x(n-1);   z2  = x(n);
        else
            th1 = yth(iDepth-1); th2 = yth(iDepth+1);
            dt1 = nth(iDepth-1); dt2 = nth(iDepth+1);
            z1  = x(iDepth-1);   z2  = x(iDepth+1);
        end
        B               = abs(z2 - z1);
        gradTh(iDepth)      = (th2 - th1) / (z2 - z1);
        deltaGradTh(iDepth) = sqrt((dt1 / B) ^ 2 + (dt2 / B) ^ 2);
    end

    %   smoothed + masked gradient from saved output ::
    statIdx      = gp15_obs.stationNo == sn & gp15_obs.depth <= maxDepth;
    gradThHat    = gp15_grad.vertGrad(statIdx);
    deltaGradThHat = gp15_grad.vertGradError(statIdx);

    %   depth markers (match doVertGradCalc station-number keying) ::
    mldVal = gp15_mld.('MLD JAK')(gp15_mld.('Station No') == sn);
    ppzVal = gp15_ppz.ppzDepth(gp15_ppz.stationNo == sn);
    eqpVal = gp15_eqp.eqp(gp15_eqp.stationNo == sn);
    zBot   = min(eqpVal, ppzVal);

    %   plot ::
    figure;

    hold('on');

    %   raw finite-difference gradient with measurement uncertainty ::
    ebRaw = errorbar(gradTh, x, deltaGradTh, deltaGradTh, 'horizontal', ...
                     'color', 'k', 'lineStyle', 'none', ...
                     'marker', 'o', 'markerSize', 10, ...
                     'markerEdgeColor', 'k', 'markerFaceColor', 'w', ...
                     'lineWidth', 1);
    ebRaw.CapSize = 0;

    %   smoothed + masked gradient with combined uncertainty ::
    ebHat = errorbar(gradThHat, x, deltaGradThHat, deltaGradThHat, 'horizontal', ...
                     'color', 'k', 'lineStyle', '-', ...
                     'marker', 'o', 'markerSize', 10, ...
                     'markerEdgeColor', 'k', 'markerFaceColor', 'k', ...
                     'lineWidth', 2);
    ebHat.CapSize = 0;

    %   depth reference lines ::
    lMld  = yline(mldVal, '--k', 'lineWidth', 1);
    lZbot = yline(zBot,   ':k',  'lineWidth', 1);

    hold('off');

    legend([ebRaw, ebHat, lMld, lZbot], ...
           {'Raw FD $\pm\,\sigma$', 'Zeroed FD $\pm\,\sigma$', 'MLD', '$z_\mathrm{bot}$'}, ...
           'interpreter', 'latex', 'fontSize', 14, 'location', 'southEast');
    box('on');
    title(['\textbf{Station ' num2str(sn) '}'], 'interpreter', 'latex', 'fontSize', 20);
    xlabel('\textbf{Gradient (dpm L$^{-1}$ m$^{-1}$)}', 'interpreter', 'latex', 'fontSize', 20);
    ylabel('\textbf{Depth (m)}', 'interpreter', 'latex', 'fontSize', 20);
    ylim([0, 200]);
    xlim([min(gradTh - deltaGradTh), max(gradTh + deltaGradTh)]);
    set(gca, 'yDir', 'reverse', 'tickLabelInterpreter', 'latex', ...
        'fontSize', 16, 'fontWeight', 'bold', 'lineWidth', 1);
    set(gcf, 'units', 'inches', 'position', [2, 2, 5, 10], ...
             'paperUnits', 'inches', 'paperSize', [5, 10]);
    exportgraphics(gcf, [plot_output_basepath 'calcVertGrad/gradPlots/grad' num2str(sn) '.pdf'], ...
                   'ContentType', 'vector');

end
close('all');

%%  end program
