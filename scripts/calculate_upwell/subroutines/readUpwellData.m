%%  read data
%   u data ::
fileName = [input_basepath 'mercator/' dataProduct '/grep_' dataProduct '_daily_uo.nc'];
variableName = ['uo_' dataProduct];
if strcmp(dataProduct, 'grep')
    fileName = [input_basepath 'mercator/' dataProduct '/grep_mean_daily_uo.nc'];
    variableName = 'uo_mean';
end
u_mercator.uo = ncread(fileName, variableName);
u_mercator.longitude = ncread(fileName, 'longitude') + 360;  % originally [-180 180] domain
u_mercator.latitude = ncread(fileName, 'latitude');                         
u_mercator.depth = ncread(fileName, 'depth');
u_mercator.time = datetime(ncread(fileName, 'time') * SEC2DAY, 'convertFrom', 'epochTime', 'epoch', '1950-01-01');  % originally days since 1950-01-01 00:00:00. time is at midnight of that day. 

%   v data ::
fileName = [input_basepath 'mercator/' dataProduct '/grep_' dataProduct '_daily_vo.nc'];
variableName = ['vo_' dataProduct];
if strcmp(dataProduct, 'grep')
    fileName = [input_basepath 'mercator/' dataProduct '/grep_mean_daily_vo.nc'];
    variableName = 'vo_mean';
end
v_mercator.vo = ncread(fileName, variableName);
v_mercator.longitude = ncread(fileName, 'longitude') + 360;  % originally [-180 180] domain
v_mercator.latitude = ncread(fileName, 'latitude');
v_mercator.depth = ncread(fileName, 'depth');
v_mercator.time = datetime(ncread(fileName, 'time') * SEC2DAY, 'convertFrom', 'epochTime', 'epoch', '1950-01-01');  % originally days since 1950-01-01 00:00:00. time is at midnight of that day. 

%   clear 
clear('fileName', 'variableName'); 

%% end subroutine
