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
	interpUpwell;  

	%   display ::
	% disp(['Done with model: ' dataProduct]); 

end

%%  compare models
%   compareModels produces two figures: the full six-product comparison and a
%   sensitivity figure with FOAM (GloSea5) and GREP excluded.
%
%   FOAM outlier — cause:
%     FOAM is a consistent outlier in equatorial upwelling relative to ECCO,
%     CGLO, GLOR, and ORAS. this is not a pipeline artefact; it reflects three
%     known properties of GloSea5/NEMO:
%       (1) NEMOVAR 3D-Var assimilates T/S/SSH but does not directly correct
%           velocities. near the equator, where geostrophic balance breaks down,
%           assimilation increments inject momentum imbalances that project onto
%           horizontal divergence and inflate the diagnosed vertical velocity.
%       (2) NEMO has documented biases in the Pacific equatorial undercurrent
%           (core too deep, too weak), directly affecting divergence and hence w.
%       (3) z*-coordinate free-surface corrections can introduce residual errors
%           in the vertically integrated w (Storto et al. 2019, Ocean Sci.).
%     GREP is excluded from the sensitivity figure because it incorporates FOAM
%     in its average and is therefore not independent of the outlier.
%
%   FOAM outlier — consequence for reported uncertainty:
%     FOAM is the primary driver of sigma_model (inter-model std across the five
%     independent products). including FOAM makes the reported total uncertainty
%     conservative — wider than the four-product agreement alone would imply.
%     this is scientifically appropriate: if FOAM's equatorial w is unreliable,
%     the correct response is a wider uncertainty band, not silent exclusion.
%
%   TODO (collation stage): when collateStations.m assembles the model input
%   table, compute and store two ensemble mean upwelling fields per station:
%     wMeanFull   — mean of ECCO, CGLO, FOAM, GLOR, ORAS (n=5, current default)
%     wMeanNoFOAM — mean of ECCO, CGLO, GLOR, ORAS       (n=4, FOAM excluded)
%   run the flux model and calcError.m with both to quantify FOAM's effect on
%   the 2D correction and total uncertainty at equatorial stations. report the
%   sensitivity in supplementary material alongside all_model_transect_noFOAM.pdf. ::
compareModels;

%%  print end
disp('Done interpolating data.')

%%  end program
