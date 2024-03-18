%%  save data 
%   .mat ::
save([pro_output_basepath 'collateData/gp15_inputs.mat'], 'gp15_inputs');

%   .xlsx ::
writetable(gp15_inputs, [pro_output_basepath 'collateData/gp15_inputs.xlsx']);

%%  end subroutine
