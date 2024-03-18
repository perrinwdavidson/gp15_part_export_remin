%%  interpData - interpolating gp15 data
%--------------------------------------------------------------------------
%%  set-up environment 
setupModel;

%%  configure
configureInterp;

%%  load data
loadInterpData;

%%  interpolate per station
stationInterp;  

%%  interpolate velocities
for dataProduct = dataProducts

	%   get product name ::
	dataProduct = dataProduct{1};

	%   interpolate ::
	interpUpwell;  % WARNING: THIS TAKES HOURS TO RUN

	%   display ::
	% disp(['Done with model: ' dataProduct]); 

end

%%  compare models
compareModels; 

%%  print end
disp('Done interpolating data.')

%%  end program
