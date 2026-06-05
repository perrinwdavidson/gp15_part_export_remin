%%  read npp
%   read in ::
%%% get directory ::
fnames = dir([input_basepath 'cbpm/']); 

%%% get length ::
n = length(fnames) - 3;  % directory information
j = 0; 

%%% preallocate ::
npp = NaN(1080, 2160, n);
nppTime = zeros(n, 1);

%%% loop and read ::
for i = 1 : 1 : length(fnames)
	if ~strcmp(fnames(i).name(1), 'c')
		continue
	else
		j = j + 1;
	end
	fname = [input_basepath 'cbpm/' fnames(i).name];
	npp(:, :, j) = hdfread(fname, 'npp', 'index', {[1  1], [1  1], [1080  2160]});
	nppTime(j) = datenum(2018, 1, 1) + str2double(extractBetween(fname, 'cbpm.2018', '.hdf')) - 1;
	%   dates of cruise are 268 to 326 or end-13 to end-6
end

%   quality control ::
NPP = double(npp);
NPP(NPP < 0) = 0;  % remove fill values 
NPP(NPP == 0) = NaN; 
% NPP = nanmean(NPP(:, :, end-13:end-6), 3) / 12;  % mg C to mmol C (averaging over cruise dates)
NPP = NPP / 12;  % mg C to mmol C

%   calculate coordinate matrices ::
[vert, hor, ~] = size(NPP);  % vert number of lat boxes, hor long boxes
latdpb = 180 / vert;  % degrees per lat box
longdpb = 360 / hor;  % degrees per long box
long_0 = -180;  % zero point
lat_0 = 90;  % zero point
lat_1 = lat_0 - (latdpb / 2);  % middle of lat box 1
long_1 = long_0 + (latdpb / 2);  % middle of long box 1
latvect = (-lat_1) : latdpb : (90 - (latdpb / 2));
lat = flipud(latvect(:));
lon = long_1 : longdpb : (180 - (longdpb / 2));
latmatrix = repmat(latvect', 1, hor); 
longmatrix = repmat(lon, vert, 1); 

%   make .nc file ::
fname = [pro_output_basepath 'readData/cbpm/npp.nc'];
if isfile(fname)
	delete(fname);
end
nccreate(fname, 'npp', 'dimensions', {'lat', length(lat), 'lon', length(lon), 'time', length(nppTime)}, 'format', 'classic');
nccreate(fname, 'lat', 'dimensions', {'lat', length(lat)}, 'format', 'classic');
nccreate(fname, 'lon', 'dimensions', {'lon', length(lon)}, 'format', 'classic');
nccreate(fname, 'time', 'dimensions', {'time', length(nppTime)}, 'format', 'classic');
ncwrite(fname, 'npp', NPP)
ncwrite(fname, 'lat', lat)
ncwrite(fname, 'lon', lon)
ncwrite(fname, 'time', nppTime)

%   save data in .mat file ::
save([pro_output_basepath 'readData/cbpm/npp.mat'], 'NPP', 'lat', 'lon', 'nppTime');

%   save data ::
writematrix(latmatrix, [pro_output_basepath 'readData/cbpm/npp_lat.csv'])
writematrix(longmatrix, [pro_output_basepath 'readData/cbpm/npp_lon.csv'])

%% end subroutine
