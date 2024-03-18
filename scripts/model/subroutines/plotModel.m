%%  plotModel - plotting model output
%%-------------------------------------------------------------------------
%%  plot 2d
figure;
tiledlayout(2, 1); 

nexttile();
    bar(gp15_flux100.stationNo, gp15_flux100.th234FluxCumul2d_ecco, 'FaceColor', 'white','EdgeColor','black','LineWidth',1.5)
    hold on 
    err = errorbar(gp15_flux100.stationNo, gp15_flux100.th234FluxCumul2d_ecco, gp15_flux100.totalTh234FluxError, 'lineWidth', 1.5); 
    err.Color = [0 0 0];                            
    err.LineStyle = 'none';
    yline(0, 'lineWidth', 1.5, 'color', 'k');
    hold off
        title('\textbf{100m}', 'interpreter', 'latex', 'fontSize', 20);
        xlabel('\textbf{Station Number}', ...
               'interpreter', 'latex', 'fontSize', 20);
        ylabel('\textbf{Flux (dpm m$^{-2}$ day$^{-1}$)}', ...
               'interpreter', 'latex', 'fontSize', 20);
        legend('Flux', 'Error', ...
               'interpreter', 'latex', 'fontSize', 20); 
            box on
            set(gca, 'tickLabelInterpreter', 'latex', 'xDir', 'reverse', ...
                'fontSize', 16, 'fontWeight', 'bold', 'lineWidth', 1);
            
nexttile();
    bar(gp15_fluxPpz.stationNo, gp15_fluxPpz.th234FluxCumul2d_ecco, 'FaceColor', 'white','EdgeColor','black','LineWidth',1.5)
    hold on 
    err = errorbar(gp15_fluxPpz.stationNo, gp15_fluxPpz.th234FluxCumul2d_ecco, gp15_fluxPpz.totalTh234FluxError, 'lineWidth', 1.5); 
    err.Color = [0 0 0];                            
    err.LineStyle = 'none';
    yline(0, 'lineWidth', 1.5, 'color', 'k');
    hold off
        title('\textbf{PPZ}', 'interpreter', 'latex', 'fontSize', 20);
        xlabel('\textbf{Station Number}', ...
               'interpreter', 'latex', 'fontSize', 20);
        ylabel('\textbf{Flux (dpm m$^{-2}$ day$^{-1}$)}', ...
               'interpreter', 'latex', 'fontSize', 20);
        legend('Flux', 'Error', ...
               'interpreter', 'latex', 'fontSize', 20); 
            box on
            set(gca, 'tickLabelInterpreter', 'latex', 'xDir', 'reverse', ...
                'fontSize', 16, 'fontWeight', 'bold', 'lineWidth', 1);
            set(gcf, 'Units', 'Inches', 'Position', [1.5, 2.5, 20, 2*6.5], ... 
                'PaperUnits', 'Inches', 'PaperSize', [20, 2*6.5])
            exportgraphics(gcf, ...
                   [plot_output_basepath 'gp15Model/gp15Plots/gp15_ecco_flux.png'], ...
                   'resolution', 300); 
            
%%  plot 2d correction
figure;
tiledlayout(2, 1); 
nexttile();
    bar(gp15_flux100.stationNo, gp15_flux100.th234FluxCumul2d_ecco - gp15_flux100.th234FluxCumul1d, ...
        'FaceColor', 'white','EdgeColor','black','LineWidth',1.5)
    hold on 
    err = errorbar(gp15_flux100.stationNo, gp15_flux100.th234FluxCumul2d_ecco - gp15_flux100.th234FluxCumul1d, ...
                   gp15_flux100.totalTh234FluxError - gp15_flux100.uncert1d, ...
                   'lineWidth', 1.5); 
    err.Color = [0 0 0];                            
    err.LineStyle = 'none';
    yline(0, 'lineWidth', 1.5, 'color', 'k');
    hold off
        title('\textbf{100m}', 'interpreter', 'latex', 'fontSize', 20);
        xlabel('\textbf{Station Number}', ...
               'interpreter', 'latex', 'fontSize', 20);
        ylabel('\textbf{Flux Correction (dpm m$^{-2}$ day$^{-1}$)}', ...
               'interpreter', 'latex', 'fontSize', 20);
        legend('Flux', 'Error', ...
               'interpreter', 'latex', 'fontSize', 20); 
            box on
            set(gca, 'tickLabelInterpreter', 'latex', 'xDir', 'reverse', ...
                'fontSize', 16, 'fontWeight', 'bold', 'lineWidth', 1);  

nexttile();
    bar(gp15_fluxPpz.stationNo, gp15_fluxPpz.th234FluxCumul2d_ecco - gp15_fluxPpz.th234FluxCumul1d, ...
        'FaceColor', 'white','EdgeColor','black','LineWidth',1.5)
    hold on 
    err = errorbar(gp15_fluxPpz.stationNo, gp15_fluxPpz.th234FluxCumul2d_ecco - gp15_fluxPpz.th234FluxCumul1d, ...
                   gp15_fluxPpz.totalTh234FluxError - gp15_fluxPpz.uncert1d, ...
                   'lineWidth', 1.5); 
    err.Color = [0 0 0];                            
    err.LineStyle = 'none';
    yline(0, 'lineWidth', 1.5, 'color', 'k');
    hold off
        title('\textbf{PPZ}', 'interpreter', 'latex', 'fontSize', 20);
        xlabel('\textbf{Station Number}', ...
               'interpreter', 'latex', 'fontSize', 20);
        ylabel('\textbf{Flux Correction (dpm m$^{-2}$ day$^{-1}$)}', ...
               'interpreter', 'latex', 'fontSize', 20);
        legend('Flux', 'Error', ...
               'interpreter', 'latex', 'fontSize', 20); 
            box on
            set(gca, 'tickLabelInterpreter', 'latex', 'xDir', 'reverse', ...
                'fontSize', 16, 'fontWeight', 'bold', 'lineWidth', 1);
            set(gcf, 'Units', 'Inches', 'Position', [1.5, 2.5, 20, 2*6.5], ... 
                'PaperUnits', 'Inches', 'PaperSize', [20, 2*6.5])
            exportgraphics(gcf, ...
                   [plot_output_basepath 'gp15Model/gp15Plots/gp15_ecco_upwellCorrect.png'], ...
                   'resolution', 300); 
            
%%  plot sensitivities
figure;
tiledlayout(2, 1); 
    nexttile();
    hold on 
    errorbar(gp15_flux100.stationNo, gp15_flux100.th234FluxCumul2d_ecco, gp15_flux100.uncertUpwell_ecco); 
    errorbar(gp15_flux100.stationNo, gp15_flux100.th234FluxCumul2d_cglo, gp15_flux100.uncertUpwell_cglo); 
    errorbar(gp15_flux100.stationNo, gp15_flux100.th234FluxCumul2d_foam, gp15_flux100.uncertUpwell_foam); 
    errorbar(gp15_flux100.stationNo, gp15_flux100.th234FluxCumul2d_glor, gp15_flux100.uncertUpwell_glor); 
    errorbar(gp15_flux100.stationNo, gp15_flux100.th234FluxCumul2d_oras, gp15_flux100.uncertUpwell_oras); 
    hold off 
        title('\textbf{100m}', 'interpreter', 'latex', 'fontSize', 20);
        xlabel('\textbf{Station Number}', ...
               'interpreter', 'latex', 'fontSize', 20);
        ylabel('\textbf{Flux (dpm m$^{-2}$ day$^{-1}$)}', ...
               'interpreter', 'latex', 'fontSize', 20); 
            set(gca, 'tickLabelInterpreter', 'latex', 'xDir', 'reverse', 'box', 'on', ...
                'fontSize', 16, 'fontWeight', 'bold', 'lineWidth', 1);
            
nexttile();
    hold on 
    errorbar(gp15_fluxPpz.stationNo, gp15_fluxPpz.th234FluxCumul2d_ecco, gp15_fluxPpz.uncertUpwell_ecco); 
    errorbar(gp15_fluxPpz.stationNo, gp15_fluxPpz.th234FluxCumul2d_cglo, gp15_fluxPpz.uncertUpwell_cglo); 
    errorbar(gp15_fluxPpz.stationNo, gp15_fluxPpz.th234FluxCumul2d_foam, gp15_fluxPpz.uncertUpwell_foam); 
    errorbar(gp15_fluxPpz.stationNo, gp15_fluxPpz.th234FluxCumul2d_glor, gp15_fluxPpz.uncertUpwell_glor); 
    errorbar(gp15_fluxPpz.stationNo, gp15_fluxPpz.th234FluxCumul2d_oras, gp15_fluxPpz.uncertUpwell_oras); 
    hold off 
        title('\textbf{PPZ}', 'interpreter', 'latex', 'fontSize', 20);
        xlabel('\textbf{Station Number}', ...
               'interpreter', 'latex', 'fontSize', 20);
        ylabel('\textbf{Flux (dpm m$^{-2}$ day$^{-1}$)}', ...
               'interpreter', 'latex', 'fontSize', 20);
        lg = legend('ECCO', 'CGLO', 'FOAM', 'GLOR', 'ORAS', ...
               'interpreter', 'latex', 'fontSize', 20); 
        lg.Layout.Tile = 'east'; 
            set(gca, 'tickLabelInterpreter', 'latex', 'xDir', 'reverse', 'box', 'on', ...
                'fontSize', 16, 'fontWeight', 'bold', 'lineWidth', 1);
            set(gcf, 'Units', 'Inches', 'Position', [1.5, 2.5, 20, 2*6.5], ... 
                'PaperUnits', 'Inches', 'PaperSize', [20, 2*6.5])
            exportgraphics(gcf, ...
                   [plot_output_basepath 'gp15Model/gp15Plots/gp15_sensitivity.png'], ...
                   'resolution', 300); 
               
%%  end routine
