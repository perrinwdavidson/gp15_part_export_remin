%%  save data
%   gp15 sampling (we qc above) ::
writetable(gp15_npp, [pro_output_basepath 'doQc/gp15/gp15_npp.xlsx'], 'writeMode', 'overwritesheet');
save([pro_output_basepath 'doQc/gp15/gp15_npp.mat'], 'gp15_npp');

%   mean npp ::
writematrix(npp_mean, [pro_output_basepath 'calcNpp/npp_mean.csv']);

%   standard error npp ::
writematrix(npp_se, [pro_output_basepath 'calcNpp/npp_se.csv']);
save([pro_output_basepath 'calcNpp/npp_se.mat'], 'npp_se');
% save([pro_output_basepath 'calcNpp/npp_se.mat'], 'npp_se');  % duplicate

%%  end subroutine
