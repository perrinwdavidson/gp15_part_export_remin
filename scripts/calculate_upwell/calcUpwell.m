%%  calcUpwell - calculating upwelling in a suite of MERCATOR models
%--------------------------------------------------------------------------
%%  set-up environment 
setupModel;

%%  load data
loadUpwellData;

%%  configure
configureUpwell;

%%  loop through all model and calculate upwelling
for dataProduct = dataProducts

	%   get product name ::
	dataProduct = dataProduct{1};

	%   read data ::
	readUpwellData;

	%   calculate upwelling ::
	calculateUpwelling;	

	%   spatially and temporally average ::
	spatialTimeAverage;

	%   plot ::
	plotUpwell;

	%   save ::
	saveUpwell;

	%   display :;
	% disp(['Done with model: ' dataProduct]); 

end

%%  print end
disp('Done calculating upwelling.')

%%  end program
