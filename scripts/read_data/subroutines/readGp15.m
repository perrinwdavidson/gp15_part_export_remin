%%  read gp15 data
%   station and cast numbers ::
fileName = [input_basepath 'gp15/gp15_stat_cast.xlsx'];
gp15_stat_cast = readtable(fileName, 'variableNamingRule', 'preserve', 'sheet', 'Sheet 1');
                          
%   observations ::
%%% total radionuclide data ::
fileName = [input_basepath 'gp15/gp15_master.xlsx'];
gp15_obs_raw = readtable(fileName, 'variableNamingRule', 'preserve', 'sheet', 'tTh-234_Summary', 'range', 'A6:AE705', 'readVariableNames', true);

%%% make array ::
gp15_obs = gp15_obs_raw(:, [2 3 6 9:12 14:17]); 

%%% give serial numbers ::
endArray = size(gp15_obs, 1); 
gp15_obs.serialNo = [1 : 1 : endArray]'; 
gp15_obs = gp15_obs(:, [end, 1:(end-1)]);

%%% particulate data ::
gp15_part_raw = readtable(fileName, 'variableNamingRule', 'preserve', 'sheet', 'pTh234', 'range', 'A2:CP439'); 

%%% make array ::
gp15_part = gp15_part_raw(:, {'stationNo', 'castNo', 'latitude', 'longitude', 'date', 'depth', 'th234PartSsf', 'uncertTh234PartSsf', 'th234PartLsf', 'uncertTh234PartLsf'});  

%   ctd ::
fileName = [input_basepath 'gp15/gp15_ctd.xlsx'];
gp15_ctd = readtable(fileName, 'variableNamingRule', 'preserve'); 

%   bottles ::
fileName = [input_basepath 'gp15/gp15_bottles.xlsx'];
gp15_bottles = readtable(fileName,  'variableNamingRule', 'preserve', 'range', 'A2:CA5126', 'readVariableNames', true); 

%   hplc data ::
fileName = [input_basepath 'gp15/gp15_hplc.xlsx'];
gp15_hplc = readtable(fileName, 'variableNamingRule', 'preserve');

%   particle data ::
fileName = [input_basepath 'gp15/gp15_particles.xlsx'];
gp15_particles = readtable(fileName, 'variableNamingRule', 'preserve');

%%  save 
%   .mat files ::
save([pro_output_basepath 'readData/gp15/gp15_obs.mat'], 'gp15_obs');
save([pro_output_basepath 'readData/gp15/gp15_part.mat'], 'gp15_part');
save([pro_output_basepath 'readData/gp15/gp15_stat_cast.mat'], 'gp15_stat_cast');
save([pro_output_basepath 'readData/gp15/gp15_bottles.mat'], 'gp15_bottles');
save([pro_output_basepath 'readData/gp15/gp15_ctd.mat'], 'gp15_ctd');
save([pro_output_basepath 'readData/gp15/gp15_hplc.mat'], 'gp15_hplc');
save([pro_output_basepath 'readData/gp15/gp15_particles.mat'], 'gp15_particles');

%   .xlsx files ::
%   n.b.: only re-write the quality controlled one (so far).
writetable(gp15_obs, [pro_output_basepath 'readData/gp15/gp15_obs.xlsx'], 'writeMode', 'overwritesheet');
writetable(gp15_part, [pro_output_basepath 'readData/gp15/gp15_part.xlsx'], 'writeMode', 'overwritesheet');

%% end subroutine
