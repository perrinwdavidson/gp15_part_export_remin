%% save output
%  excel ::
mldCalcsFinal = NaN(NUMSTAT, 6); 
mldCalcsFinal(:, 1) = gp15_stations.stationNo'; 
mldCalcsFinal(:, 2) = mld_JAK;
mldCalcsFinal(:, 3) = mld_dBM;
mldCalcsFinal(:, 4) = mld_BW(:, 1);
mldCalcsFinal(:, 5) = mld_BW(:, 2);
mldCalcsFinal(:, 6) = mld_BW(:, 3);
mldCalcsFinal = array2table(mldCalcsFinal);
mldCalcsFinal.Properties.VariableNames = {'Station No', 'MLD JAK', 'MLD dBM', 'MLD BW 01', ...
                                          'MLD BW 05', 'MLD BW 125'};
writetable(mldCalcsFinal, [sim_output_basepath 'calcMld/gp15_mld.xlsx'])

%   mat ::
gp15_mld = mldCalcsFinal;
save([sim_output_basepath 'calcMld/gp15_mld.mat'], 'gp15_mld');

%%  clean up
clear('mldCalcsFinal'); 
close('all'); 

%% end subroutine
