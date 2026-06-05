%%  collectData - collecting and collating all model data inputs
%--------------------------------------------------------------------------
%%  set-up environment 
setupModel;

%%  load data
loadCollateData;

%%  collect data per stations
%   collateStations.m joins per-product kriged upwelling (w_ecco, w_cglo,
%   w_foam, w_glor, w_oras, w_grep) as individual columns in gp15_inputs,
%   then computes two ensemble mean columns and their propagated kriging
%   uncertainties:
%
%     w_meanFull / wErr_meanFull     — mean of 5 independent products
%                                      (ECCO, CGLO, FOAM, GLOR, ORAS)
%                                      wErr = sqrt(Σ wErr_k²) / 5
%     w_meanNoFOAM / wErr_meanNoFOAM — mean of 4 products excluding FOAM
%                                      (ECCO, CGLO, GLOR, ORAS)
%                                      wErr = sqrt(Σ wErr_k²) / 4
%
%   calcError.m computes totalXxxFluxError_meanFull and _meanNoFOAM by
%   combining sigma_2d (from the respective mean w) with sigma_model (inter-
%   model std of 5 or 4 members respectively) in quadrature. sigma_model is
%   attached only to the mean-ensemble products; individual members carry
%   sigma_2d only.
%
%   calcStats.m exposes fluxProduct / uncertProduct variables (default both
%   'meanFull') to select which flux and which total error are reported in
%   depthData and the manuscript tables. ::
collateStations;

%%  save data 
saveCollateData;

%%  print
disp('Done collating data.');

%%  end program
