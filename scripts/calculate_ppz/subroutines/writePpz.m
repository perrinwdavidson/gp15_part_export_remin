%%  write out array ::
%   .mat ::
save([sim_output_basepath 'calcPpz/gp15/gp15_ppz.mat'], 'gp15_ppz');

%   .xlsx ::
writetable(gp15_ppz, [sim_output_basepath 'calcPpz/gp15/gp15_ppz.xlsx']);

%%  end subroutine
