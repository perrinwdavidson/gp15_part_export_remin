%%  configure
%   set averaging time ::
deltaT = 16;  % set as twice the time resolution

%   set initial length bounds ::
lxCurrent = lx;
ltCurrent = deltaT; 

%   set tolerances ::
XTOL = 1E-10;
TTOL = 1E-1;

%%  end subroutine
