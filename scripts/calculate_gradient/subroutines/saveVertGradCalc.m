%%  save data
matPath = [sim_output_basepath 'calcVertGrad/gp15_grad.mat'];
save(matPath, 'gp15_grad');
writeHash(matPath);

%%  end subroutine
