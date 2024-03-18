%%  read data
load([sim_output_basepath 'calcPpz/gp15/gp15_ppz.mat'], 'gp15_ppz');

%%  sort rows
gp15_ppz = sortrows(gp15_ppz, 'stationNo');

%%  remove misc and non-integer stations
gp15_ppz = rmmissing(gp15_ppz);
gp15_ppz = gp15_ppz(mod(gp15_ppz.stationNo, 1) == 0, :);

%%  save data
%   excel ::
writetable(gp15_ppz, [pro_output_basepath 'doQc/ppz/gp15_ppz.xlsx']);

%   mat ::
save([pro_output_basepath 'doQc/ppz/gp15_ppz.mat'], 'gp15_ppz');

%% end subroutine
