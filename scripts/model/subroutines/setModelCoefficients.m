%%  set coefficients
TH_HALFLIFE = 24.101;
LAMBDA = log(2.000) / TH_HALFLIFE;

%%  set spline smoothing parameters
%   p = 0: maximally smooth (straight line); p = 1: interpolating spline.
%   run sensitivity at p_nominal, p/2, and p*2 before finalising values. ::
P_SPLINE_GRADIENT = 0.9;   % sparse/noisy profiles: Th-234 gradient, Th234/U238 ratio
P_SPLINE_DENSITY  = 0.99;  % dense/precise CTD profiles: potential density, temperature

%%  end subroutine
