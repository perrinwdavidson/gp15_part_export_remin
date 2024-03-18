%%  configure
%   at what depth do you want the profile calculations and plotting to end?
BTM_DEPTH = 400;

%   how far below mld guess do you want to plot?
ADD_DEPTH = 50;

%   difference for potential density anamoly?
MAKE_ANAMOLY = 1000;

%   what dBM reference depth and threshold -- from dBM (2004, AGU:O)?
DELTA_DEPTH = 10;  

%   what is the delta temperature -- from dBM (2004, AGU:O)?
DELTA_T = 0.2;

%   what thresholds for BW method -- from BW (2009, GBC)? *
DELTA_POTDENS = [0.01, 0.05, 0.125]; 

%   what averaging depth for surface BW value? 
MLD_AVE_DEPTH = 10; 

%   what standard deviation level for Pickart method -- from Pickart (2002, JPO)?
NUM_STD = 2;

%   do you want to plott?
plotting = 'yes';

%--------------------------------------------------------------------------
% * it is most clearly written in the 12th reference/notes of Bishop et al. 
%   (2004, Science)
%%  end subroutine
