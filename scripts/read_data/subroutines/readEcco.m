%%  read ecco upwelling data
%   get filename ::
fileName = [input_basepath 'ecco/ecco_wo_quart1.cdf'];

%   coordinates ::
ecco_wo.latitude = ncread(fileName, 'lat');
ecco_wo.longitude = ncread(fileName, 'lon');
ecco_wo.depth = ncread(fileName, 'depth');

%  get time and data ::
iStart = 1;
for i = iStart : 1 : 4

	% get filename ::
	fileName = [input_basepath 'ecco/ecco_wo_quart' num2str(i) '.cdf'];

	% get time ::
	if i == iStart
		ecco_wo.time = ncread(fileName, 'time'); 
		ecco_wo.wo = ncread(fileName, 'Wave');
	else
		ecco_wo.time = cat(1, ecco_wo.time, ncread(fileName, 'time')); 
		ecco_wo.wo = cat(4, ecco_wo.wo, ncread(fileName, 'Wave'));
	end
end

%%  quality control
ecco_wo.wo(ecco_wo.wo == -1.000000000000000e+10) = NaN;
ecco_wo.time = ecco_wo.time .* (60 ^ 2);
ecco_wo.time = datetime(ecco_wo.time, 'convertFrom', 'posixtime');

%%  save
%   n.b.: do not save to .nc as will not need to review (literally just what is above).
save([pro_output_basepath 'readData/ecco/ecco_wo.mat'], 'ecco_wo');

%% end subroutine
