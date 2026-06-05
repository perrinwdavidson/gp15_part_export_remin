%%  save data
%   .mat ::
matPath = [pro_output_basepath 'collateData/gp15_inputs.mat'];
save(matPath, 'gp15_inputs');
writeHash(matPath);

%   .xlsx ::
writetable(gp15_inputs, [pro_output_basepath 'collateData/gp15_inputs.xlsx']);

%%  end subroutine
