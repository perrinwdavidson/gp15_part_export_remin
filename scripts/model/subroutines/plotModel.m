%%  plotModel - plotting model output
%%------------------------------------------------------------------------

%%  shared settings
products   = {'ecco', 'cglo', 'foam', 'glor', 'oras', 'grep', 'meanFull', 'meanNoFOAM'};
prodLabels = {'ECCO', 'CGLO', 'FOAM', 'GLOR', 'ORAS', 'GREP', 'Mean (Full)', 'Mean (No FOAM)'};

depthTbls   = {gp15_flux100, gp15_fluxPpz};
depthTags   = {'100m', 'ppz'};
depthTitles = {'100 m', 'PPZ'};

rowLabels    = {'1D Flux', 'Upwelling Correction', '2D Flux'};
th234Ylabels = {'$^{234}$Th Flux (dpm m$^{-2}$ d$^{-1}$)', ...
                '$^{234}$Th Correction (dpm m$^{-2}$ d$^{-1}$)', ...
                '$^{234}$Th Flux (dpm m$^{-2}$ d$^{-1}$)'};
pocYlabels   = {'POC Flux (mmol m$^{-2}$ d$^{-1}$)', ...
                'POC Correction (mmol m$^{-2}$ d$^{-1}$)', ...
                'POC Flux (mmol m$^{-2}$ d$^{-1}$)'};

eccoColor   = [0.00, 0.45, 0.74];
cmemsColor  = [0.65, 0.65, 0.65];
cmemsStyles = {'-', '--', '-.', ':'};
meanColor   = [0, 0, 0];
grepColor   = [0.85, 0.10, 0.10];

figW = 20;  figH = 18;
outDir = [plot_output_basepath 'gp15Model/gp15Plots/'];

%%  1. per-product bar plots (3 rows x 2 cols; one figure per product per depth)
for iProd = 1 : length(products)
    prod = products{iProd};
    pLab = prodLabels{iProd};

    for iDepth = 1 : 2
        tbl  = depthTbls{iDepth};
        dTag = depthTags{iDepth};
        dLab = depthTitles{iDepth};
        stn  = tbl.stationNo;

        th1d  = tbl.th234FluxCumul1d;
        poc1d = tbl.pocFluxCumul1d;

        thData  = {th1d, ...
                   tbl.(['th234FluxCumul2d_' prod]) - th1d, ...
                   tbl.(['th234FluxCumul2d_' prod])};
        pocData = {poc1d, ...
                   tbl.(['pocFluxCumul2d_' prod]) - poc1d, ...
                   tbl.(['pocFluxCumul2d_' prod])};
        thErr   = {tbl.uncert1d, ...
                   abs(tbl.(['uncertUpwellCorrect_'    prod])), ...
                   tbl.(['uncertUpwell_'               prod])};
        pocErr  = {tbl.uncert1dPoc, ...
                   abs(tbl.(['uncertUpwellPocCorrect_' prod])), ...
                   tbl.(['uncertUpwellPoc_'            prod])};

        figure;
        tl = tiledlayout(3, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
        title(tl, ['\textbf{' pLab ' -- ' dLab '}'], ...
              'interpreter', 'latex', 'fontSize', 20);

        for iRow = 1 : 3
            for iCol = 1 : 2
                nexttile((iRow - 1) * 2 + iCol);
                if iCol == 1
                    yv = thData{iRow};  ev = thErr{iRow};  yl = th234Ylabels{iRow};
                else
                    yv = pocData{iRow};  ev = pocErr{iRow};  yl = pocYlabels{iRow};
                end
                bar(stn, yv, 'FaceColor', [1 1 1], 'EdgeColor', [0 0 0], 'LineWidth', 1.5);
                hold('on');
                errorbar(stn, yv, ev, 'Color', [0 0 0], 'LineStyle', 'none', ...
                         'LineWidth', 1.5, 'CapSize', 0);
                yline(0, 'Color', 'k', 'LineWidth', 1.5);
                hold('off');
                title(['\textbf{' rowLabels{iRow} '}'], 'interpreter', 'latex', 'fontSize', 16);
                ylabel(['\textbf{' yl '}'], 'interpreter', 'latex', 'fontSize', 13);
                if iRow == 3
                    xlabel('\textbf{Station Number}', 'interpreter', 'latex', 'fontSize', 13);
                end
                box('on');
                set(gca, 'tickLabelInterpreter', 'latex', 'xDir', 'reverse', ...
                    'fontSize', 12, 'fontWeight', 'bold', 'lineWidth', 1);
            end
        end

        set(gcf, 'Units', 'Inches', 'Position', [1, 1, figW, figH], ...
            'PaperUnits', 'Inches', 'PaperSize', [figW, figH]);
        exportgraphics(gcf, [outDir 'gp15_' prod '_' dTag '.pdf'], 'ContentType', 'vector');
    end
end

%%  2-4. comparison line plots (one figure per ensemble type per depth)
%   POC flux is NaN at stations without valid POC:Th ratio data (no particulate
%   sampling or all observations flagged by QC). Filter per-line with ~isnan so
%   each product's line is continuous over the stations where it has data.
for iDepth = 1 : 2
    tbl  = depthTbls{iDepth};
    dTag = depthTags{iDepth};
    dLab = depthTitles{iDepth};
    stn  = tbl.stationNo;

    th1d  = tbl.th234FluxCumul1d;
    poc1d = tbl.pocFluxCumul1d;
    e1d   = tbl.uncert1d;
    e1dP  = tbl.uncert1dPoc;

    %%  2. meanFull ensemble: ECCO (blue) + CGLO/FOAM/GLOR/ORAS (gray) + mean (thick black)
    memFull     = {'ecco', 'cglo', 'foam', 'glor', 'oras'};
    memFullLabs = {'ECCO', 'CGLO', 'FOAM', 'GLOR', 'ORAS'};

    figure;
    tl = tiledlayout(3, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tl, ['\textbf{Full Ensemble Members vs Mean -- ' dLab '}'], ...
          'interpreter', 'latex', 'fontSize', 20);

    for iRow = 1 : 3
        for iCol = 1 : 2
            nexttile((iRow - 1) * 2 + iCol);
            hold('on');
            if iCol == 1;  yl = th234Ylabels{iRow};  else;  yl = pocYlabels{iRow};  end

            if iRow == 1
                if iCol == 1;  yv = th1d;  ev = e1d;  else;  yv = poc1d;  ev = e1dP;  end
                ok = ~isnan(yv) & ~isnan(ev);
                eb = errorbar(stn(ok), yv(ok), ev(ok), '-', 'Color', meanColor, 'LineWidth', 2, 'CapSize', 0);
                eb.DisplayName = '1D';
            else
                for iMem = 1 : 5
                    m = memFull{iMem};
                    if iCol == 1
                        if iRow == 2
                            yv = tbl.(['th234FluxCumul2d_' m]) - th1d;
                            ev = abs(tbl.(['uncertUpwellCorrect_' m]));
                        else
                            yv = tbl.(['th234FluxCumul2d_' m]);
                            ev = tbl.(['uncertUpwell_' m]);
                        end
                    else
                        if iRow == 2
                            yv = tbl.(['pocFluxCumul2d_' m]) - poc1d;
                            ev = abs(tbl.(['uncertUpwellPocCorrect_' m]));
                        else
                            yv = tbl.(['pocFluxCumul2d_' m]);
                            ev = tbl.(['uncertUpwellPoc_' m]);
                        end
                    end
                    if strcmp(m, 'ecco')
                        lc = eccoColor;  ls = '-';  lw = 1;
                    else
                        lc = cmemsColor;  ls = cmemsStyles{iMem - 1};  lw = 1;
                    end
                    ok = ~isnan(yv) & ~isnan(ev);
                    eb = errorbar(stn(ok), yv(ok), ev(ok), ls, 'Color', lc, 'LineWidth', lw, 'CapSize', 0);
                    eb.DisplayName = memFullLabs{iMem};
                end
                % ensemble mean
                if iCol == 1
                    if iRow == 2
                        ym = tbl.th234FluxCumul2d_meanFull - th1d;
                        em = abs(tbl.uncertUpwellCorrect_meanFull);
                    else
                        ym = tbl.th234FluxCumul2d_meanFull;
                        em = tbl.uncertUpwell_meanFull;
                    end
                else
                    if iRow == 2
                        ym = tbl.pocFluxCumul2d_meanFull - poc1d;
                        em = abs(tbl.uncertUpwellPocCorrect_meanFull);
                    else
                        ym = tbl.pocFluxCumul2d_meanFull;
                        em = tbl.uncertUpwellPoc_meanFull;
                    end
                end
                ok = ~isnan(ym) & ~isnan(em);
                eb = errorbar(stn(ok), ym(ok), em(ok), '-', 'Color', meanColor, 'LineWidth', 2.5, 'CapSize', 0);
                eb.DisplayName = 'Mean (Full)';
            end

            yline(0, 'Color', 'k', 'LineWidth', 1, 'HandleVisibility', 'off');
            hold('off');
            title(['\textbf{' rowLabels{iRow} '}'], 'interpreter', 'latex', 'fontSize', 16);
            ylabel(['\textbf{' yl '}'], 'interpreter', 'latex', 'fontSize', 13);
            if iRow == 3
                xlabel('\textbf{Station Number}', 'interpreter', 'latex', 'fontSize', 13);
            end
            if iRow == 3 && iCol == 2
                legend('location', 'best', 'interpreter', 'latex', 'fontSize', 12);
            end
            box('on');
            set(gca, 'tickLabelInterpreter', 'latex', 'xDir', 'reverse', ...
                'fontSize', 12, 'fontWeight', 'bold', 'lineWidth', 1);
        end
    end

    set(gcf, 'Units', 'Inches', 'Position', [1, 1, figW, figH], ...
        'PaperUnits', 'Inches', 'PaperSize', [figW, figH]);
    exportgraphics(gcf, [outDir 'gp15_comparison_meanFull_' dTag '.pdf'], 'ContentType', 'vector');

    %%  3. meanNoFOAM ensemble: ECCO (blue) + CGLO/GLOR/ORAS (gray) + mean (thick black)
    memNoFoam     = {'ecco', 'cglo', 'glor', 'oras'};
    memNoFoamLabs = {'ECCO', 'CGLO', 'GLOR', 'ORAS'};

    figure;
    tl = tiledlayout(3, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tl, ['\textbf{No-FOAM Ensemble Members vs Mean -- ' dLab '}'], ...
          'interpreter', 'latex', 'fontSize', 20);

    for iRow = 1 : 3
        for iCol = 1 : 2
            nexttile((iRow - 1) * 2 + iCol);
            hold('on');
            if iCol == 1;  yl = th234Ylabels{iRow};  else;  yl = pocYlabels{iRow};  end

            if iRow == 1
                if iCol == 1;  yv = th1d;  ev = e1d;  else;  yv = poc1d;  ev = e1dP;  end
                ok = ~isnan(yv) & ~isnan(ev);
                eb = errorbar(stn(ok), yv(ok), ev(ok), '-', 'Color', meanColor, 'LineWidth', 2, 'CapSize', 0);
                eb.DisplayName = '1D';
            else
                cIdx = 0;
                for iMem = 1 : 4
                    m = memNoFoam{iMem};
                    if iCol == 1
                        if iRow == 2
                            yv = tbl.(['th234FluxCumul2d_' m]) - th1d;
                            ev = abs(tbl.(['uncertUpwellCorrect_' m]));
                        else
                            yv = tbl.(['th234FluxCumul2d_' m]);
                            ev = tbl.(['uncertUpwell_' m]);
                        end
                    else
                        if iRow == 2
                            yv = tbl.(['pocFluxCumul2d_' m]) - poc1d;
                            ev = abs(tbl.(['uncertUpwellPocCorrect_' m]));
                        else
                            yv = tbl.(['pocFluxCumul2d_' m]);
                            ev = tbl.(['uncertUpwellPoc_' m]);
                        end
                    end
                    if strcmp(m, 'ecco')
                        lc = eccoColor;  ls = '-';  lw = 1;
                    else
                        cIdx = cIdx + 1;
                        lc = cmemsColor;  ls = cmemsStyles{cIdx};  lw = 1;
                    end
                    ok = ~isnan(yv) & ~isnan(ev);
                    eb = errorbar(stn(ok), yv(ok), ev(ok), ls, 'Color', lc, 'LineWidth', lw, 'CapSize', 0);
                    eb.DisplayName = memNoFoamLabs{iMem};
                end
                % ensemble mean
                if iCol == 1
                    if iRow == 2
                        ym = tbl.th234FluxCumul2d_meanNoFOAM - th1d;
                        em = abs(tbl.uncertUpwellCorrect_meanNoFOAM);
                    else
                        ym = tbl.th234FluxCumul2d_meanNoFOAM;
                        em = tbl.uncertUpwell_meanNoFOAM;
                    end
                else
                    if iRow == 2
                        ym = tbl.pocFluxCumul2d_meanNoFOAM - poc1d;
                        em = abs(tbl.uncertUpwellPocCorrect_meanNoFOAM);
                    else
                        ym = tbl.pocFluxCumul2d_meanNoFOAM;
                        em = tbl.uncertUpwellPoc_meanNoFOAM;
                    end
                end
                ok = ~isnan(ym) & ~isnan(em);
                eb = errorbar(stn(ok), ym(ok), em(ok), '-', 'Color', meanColor, 'LineWidth', 2.5, 'CapSize', 0);
                eb.DisplayName = 'Mean (No FOAM)';
            end

            yline(0, 'Color', 'k', 'LineWidth', 1, 'HandleVisibility', 'off');
            hold('off');
            title(['\textbf{' rowLabels{iRow} '}'], 'interpreter', 'latex', 'fontSize', 16);
            ylabel(['\textbf{' yl '}'], 'interpreter', 'latex', 'fontSize', 13);
            if iRow == 3
                xlabel('\textbf{Station Number}', 'interpreter', 'latex', 'fontSize', 13);
            end
            if iRow == 3 && iCol == 2
                legend('location', 'best', 'interpreter', 'latex', 'fontSize', 12);
            end
            box('on');
            set(gca, 'tickLabelInterpreter', 'latex', 'xDir', 'reverse', ...
                'fontSize', 12, 'fontWeight', 'bold', 'lineWidth', 1);
        end
    end

    set(gcf, 'Units', 'Inches', 'Position', [1, 1, figW, figH], ...
        'PaperUnits', 'Inches', 'PaperSize', [figW, figH]);
    exportgraphics(gcf, [outDir 'gp15_comparison_meanNoFOAM_' dTag '.pdf'], 'ContentType', 'vector');

    %%  4. CMEMS members, CMEMS mean, and GREP
    %   CMEMS mean computed inline as (CGLO+FOAM+GLOR+ORAS)/4;
    %   uncertainty propagated as sqrt(sum sigma_k^2)/4.
    cmemsMembers = {'cglo', 'foam', 'glor', 'oras'};
    cmemsLabs    = {'CGLO', 'FOAM', 'GLOR', 'ORAS'};

    figure;
    tl = tiledlayout(3, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tl, ['\textbf{CMEMS Members, Mean, and GREP -- ' dLab '}'], ...
          'interpreter', 'latex', 'fontSize', 20);

    for iRow = 1 : 3
        for iCol = 1 : 2
            nexttile((iRow - 1) * 2 + iCol);
            hold('on');
            if iCol == 1;  yl = th234Ylabels{iRow};  else;  yl = pocYlabels{iRow};  end

            if iRow == 1
                if iCol == 1;  yv = th1d;  ev = e1d;  else;  yv = poc1d;  ev = e1dP;  end
                ok = ~isnan(yv) & ~isnan(ev);
                eb = errorbar(stn(ok), yv(ok), ev(ok), '-', 'Color', meanColor, 'LineWidth', 2, 'CapSize', 0);
                eb.DisplayName = '1D';
            else
                % CMEMS members
                for iMem = 1 : 4
                    m = cmemsMembers{iMem};
                    if iCol == 1
                        if iRow == 2
                            yv = tbl.(['th234FluxCumul2d_' m]) - th1d;
                            ev = abs(tbl.(['uncertUpwellCorrect_' m]));
                        else
                            yv = tbl.(['th234FluxCumul2d_' m]);
                            ev = tbl.(['uncertUpwell_' m]);
                        end
                    else
                        if iRow == 2
                            yv = tbl.(['pocFluxCumul2d_' m]) - poc1d;
                            ev = abs(tbl.(['uncertUpwellPocCorrect_' m]));
                        else
                            yv = tbl.(['pocFluxCumul2d_' m]);
                            ev = tbl.(['uncertUpwellPoc_' m]);
                        end
                    end
                    ok = ~isnan(yv) & ~isnan(ev);
                    eb = errorbar(stn(ok), yv(ok), ev(ok), cmemsStyles{iMem}, ...
                                  'Color', cmemsColor, 'LineWidth', 1, 'CapSize', 0);
                    eb.DisplayName = cmemsLabs{iMem};
                end
                % CMEMS mean (CGLO+FOAM+GLOR+ORAS)/4
                if iCol == 1
                    if iRow == 2
                        ym = (tbl.th234FluxCumul2d_cglo + tbl.th234FluxCumul2d_foam + ...
                              tbl.th234FluxCumul2d_glor + tbl.th234FluxCumul2d_oras) / 4 - th1d;
                        em = sqrt(tbl.uncertUpwellCorrect_cglo .^ 2 + tbl.uncertUpwellCorrect_foam .^ 2 + ...
                                  tbl.uncertUpwellCorrect_glor .^ 2 + tbl.uncertUpwellCorrect_oras .^ 2) / 4;
                    else
                        ym = (tbl.th234FluxCumul2d_cglo + tbl.th234FluxCumul2d_foam + ...
                              tbl.th234FluxCumul2d_glor + tbl.th234FluxCumul2d_oras) / 4;
                        em = sqrt(tbl.uncertUpwell_cglo .^ 2 + tbl.uncertUpwell_foam .^ 2 + ...
                                  tbl.uncertUpwell_glor .^ 2 + tbl.uncertUpwell_oras .^ 2) / 4;
                    end
                else
                    if iRow == 2
                        ym = (tbl.pocFluxCumul2d_cglo + tbl.pocFluxCumul2d_foam + ...
                              tbl.pocFluxCumul2d_glor + tbl.pocFluxCumul2d_oras) / 4 - poc1d;
                        em = sqrt(tbl.uncertUpwellPocCorrect_cglo .^ 2 + tbl.uncertUpwellPocCorrect_foam .^ 2 + ...
                                  tbl.uncertUpwellPocCorrect_glor .^ 2 + tbl.uncertUpwellPocCorrect_oras .^ 2) / 4;
                    else
                        ym = (tbl.pocFluxCumul2d_cglo + tbl.pocFluxCumul2d_foam + ...
                              tbl.pocFluxCumul2d_glor + tbl.pocFluxCumul2d_oras) / 4;
                        em = sqrt(tbl.uncertUpwellPoc_cglo .^ 2 + tbl.uncertUpwellPoc_foam .^ 2 + ...
                                  tbl.uncertUpwellPoc_glor .^ 2 + tbl.uncertUpwellPoc_oras .^ 2) / 4;
                    end
                end
                ok = ~isnan(ym) & ~isnan(em);
                eb = errorbar(stn(ok), ym(ok), em(ok), '-', 'Color', meanColor, 'LineWidth', 2.5, 'CapSize', 0);
                eb.DisplayName = 'CMEMS Mean';
                % GREP
                if iCol == 1
                    if iRow == 2
                        yg = tbl.th234FluxCumul2d_grep - th1d;
                        eg = abs(tbl.uncertUpwellCorrect_grep);
                    else
                        yg = tbl.th234FluxCumul2d_grep;
                        eg = tbl.uncertUpwell_grep;
                    end
                else
                    if iRow == 2
                        yg = tbl.pocFluxCumul2d_grep - poc1d;
                        eg = abs(tbl.uncertUpwellPocCorrect_grep);
                    else
                        yg = tbl.pocFluxCumul2d_grep;
                        eg = tbl.uncertUpwellPoc_grep;
                    end
                end
                ok = ~isnan(yg) & ~isnan(eg);
                eb = errorbar(stn(ok), yg(ok), eg(ok), '--', 'Color', grepColor, 'LineWidth', 2, 'CapSize', 0);
                eb.DisplayName = 'GREP';
            end

            yline(0, 'Color', 'k', 'LineWidth', 1, 'HandleVisibility', 'off');
            hold('off');
            title(['\textbf{' rowLabels{iRow} '}'], 'interpreter', 'latex', 'fontSize', 16);
            ylabel(['\textbf{' yl '}'], 'interpreter', 'latex', 'fontSize', 13);
            if iRow == 3
                xlabel('\textbf{Station Number}', 'interpreter', 'latex', 'fontSize', 13);
            end
            if iRow == 3 && iCol == 2
                legend('location', 'best', 'interpreter', 'latex', 'fontSize', 12);
            end
            box('on');
            set(gca, 'tickLabelInterpreter', 'latex', 'xDir', 'reverse', ...
                'fontSize', 12, 'fontWeight', 'bold', 'lineWidth', 1);
        end
    end

    set(gcf, 'Units', 'Inches', 'Position', [1, 1, figW, figH], ...
        'PaperUnits', 'Inches', 'PaperSize', [figW, figH]);
    exportgraphics(gcf, [outDir 'gp15_comparison_cmems_grep_' dTag '.pdf'], 'ContentType', 'vector');

end % iDepth

%%  end routine
