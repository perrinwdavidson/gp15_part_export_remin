%%  save data
%   excel ::
writetable(gp15_flux, [sim_output_basepath 'gp15Model/modelOutput/gp15_flux.xlsx']);
writetable(gp15_flux100, [sim_output_basepath 'gp15Model/modelOutput/gp15_flux100.xlsx']);
writetable(gp15_fluxPpz, [sim_output_basepath 'gp15Model/modelOutput/gp15_fluxPpz.xlsx']);
writetable(gp15_flux100Ppz, [sim_output_basepath 'gp15Model/modelOutput/gp15_flux100Ppz.xlsx']);

%   mat ::
fluxPaths = { ...
    [sim_output_basepath 'gp15Model/modelOutput/gp15_flux.mat'], ...
    [sim_output_basepath 'gp15Model/modelOutput/gp15_flux100.mat'], ...
    [sim_output_basepath 'gp15Model/modelOutput/gp15_fluxPpz.mat'], ...
    [sim_output_basepath 'gp15Model/modelOutput/gp15_flux100Ppz.mat'] };
save(fluxPaths{1}, 'gp15_flux');          writeHash(fluxPaths{1});
save(fluxPaths{2}, 'gp15_flux100');       writeHash(fluxPaths{2});
save(fluxPaths{3}, 'gp15_fluxPpz');       writeHash(fluxPaths{3});
save(fluxPaths{4}, 'gp15_flux100Ppz');    writeHash(fluxPaths{4});

%%  end subroutine
