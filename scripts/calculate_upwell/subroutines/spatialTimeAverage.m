%%  spatial averaging
%   specify spatial averaging ::
spatAve = floor(degreeAve / spaceResolution);

%   degrees in the longitude ::
wSpatAve = movmean(w, spatAve, 1, 'omitnan');

%   degrees in latitude ::
wSpatAve = movmean(wSpatAve, spatAve, 2, 'omitnan');

%%  temporal averaging
%   specify temporal averaging ::
timeAve = floor(deltaDay / timeResolution);

%   days in time ::
wSpatAve = movmean(w, timeAve, 4, 'omitnan');

%%  end subroutine
