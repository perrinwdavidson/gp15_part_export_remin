%%  save
%   make variables ::
wLat = u_mercator.latitude';
wLon = u_mercator.longitude';
wDepth = u_mercator.depth';
wTime = u_mercator.time';
% wFull = w;  % previous: saved raw divergence-integral; temporal averaging not applied to saved field
wFull = wSpatAve;  % spatially and temporally averaged (deltaDay = 35 d, matching Th-234 mean lifetime)

%   write ::
%%% make ncdf
fileName = [sim_output_basepath 'calcUpwell/w_' dataProduct '.nc'];
delete(fileName); 
nccreate(fileName, 'w',...
         'Dimensions', {'longitude', NUMLON, ...
                        'latitude', NUMLAT, ...
                        'depth', NUMDEPTH, ...
                        'time', NUMTIME},...
         'FillValue', NaN);  
     
nccreate(fileName, 'latitude',...
         'Dimensions', {'latitude', NUMLAT},...
         'FillValue', NaN); 
     
nccreate(fileName, 'longitude',...
         'Dimensions', {'longitude', NUMLON},...
         'FillValue', NaN);      
   
nccreate(fileName, 'depth',...
         'Dimensions', {'depth', NUMDEPTH},...
         'FillValue', NaN); 
     
nccreate(fileName, 'time',...
         'Dimensions', {'time', NUMTIME},...
         'FillValue', NaN);      
     
%%% write ncdf ::
ncwrite(fileName, 'w', wFull);
ncwrite(fileName, 'latitude', wLat);
ncwrite(fileName, 'longitude', wLon);
ncwrite(fileName, 'depth', wDepth);
ncwrite(fileName, 'time', datenum(wTime));  % n.b.: datenum

%%  end subroutine
