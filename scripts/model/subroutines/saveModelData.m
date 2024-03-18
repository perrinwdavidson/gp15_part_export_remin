%%  save data
%   excel ::
writetable(gp15_flux, [sim_output_basepath 'gp15Model/modelOutput/gp15_flux.xlsx']);
writetable(gp15_flux100, [sim_output_basepath 'gp15Model/modelOutput/gp15_flux100.xlsx']);
writetable(gp15_fluxPpz, [sim_output_basepath 'gp15Model/modelOutput/gp15_fluxPpz.xlsx']);
writetable(gp15_flux100Ppz, [sim_output_basepath 'gp15Model/modelOutput/gp15_flux100Ppz.xlsx']);

%   mat ::
save([sim_output_basepath 'gp15Model/modelOutput/gp15_flux.mat'], 'gp15_flux');
save([sim_output_basepath 'gp15Model/modelOutput/gp15_flux100.mat'], 'gp15_flux100');
save([sim_output_basepath 'gp15Model/modelOutput/gp15_fluxPpz.mat'], 'gp15_fluxPpz');
save([sim_output_basepath 'gp15Model/modelOutput/gp15_flux100Ppz.mat'], 'gp15_flux100Ppz');
save([sim_output_basepath 'gp15Model/modelOutput/textStats.mat'], 'textStats');

%%  end subroutine
