% see plotVertGrad.m for diagnostic figures

%   set max depth ::
maxDepth = 10000; 

%   make output array ::
gp15_grad = gp15_obs(:, {'stationNo', 'longitude', 'latitude', 'depth'}); 
gp15_grad.vertGrad = NaN(size(gp15_obs, 1), 1);
gp15_grad.vertGradError = NaN(size(gp15_obs, 1), 1);

%   loop through all stations ::
for iStat = 1 : 1 : NUMSTAT
    
    % find station:
    sn = gp15_stations.stationNo(iStat);

    %   get station number ::
    % statNo = gp15_stations.stationNo(iStat);  % redundant: identical to sn above
    
    % find ppz:
    ez = gp15_ppz.ppzDepth(gp15_ppz.stationNo == sn);

    % find mld ::
    % S3: use station-number keying, not positional iStat, so all three depth refs are associative ::
    % mldVal = gp15_mld.("MLD JAK")(iStat);
    mldVal = gp15_mld.("MLD JAK")(gp15_mld.("Station No") == sn);
    
    % roots ::
    rootsFt = gp15_eqp.eqp(gp15_eqp.stationNo == sn);
        
    %   get data ::
    x = gp15_obs.depth(gp15_obs.stationNo == sn &  gp15_obs.depth <= maxDepth); 
    yth = gp15_obs.th234(gp15_obs.stationNo == sn &  gp15_obs.depth <= maxDepth); 
    nth = gp15_obs.uncertTh234(gp15_obs.stationNo == sn &  gp15_obs.depth <= maxDepth); 
    
    %   calculate number of depths ::
    n = length(x);
    gradTh      = zeros(n, 1);
    deltaGradTh = zeros(n, 1);

    %   calculate gradient ::
    for iDepth = 1 : 1 : n
        
        %   calculate for forward difference ::
        if iDepth == 1
            
            %   get thorim values ::
            th1 = yth(iDepth); 
            th2 = yth(iDepth + 1); 
            
            %   get thorium error values ::
            deltaTh1 = nth(iDepth); 
            deltaTh2 = nth(iDepth + 1); 
            
            %   get depth values ::
            z1 = x(iDepth); 
            z2 = x(iDepth + 1);
            
            %   calculate B ::
            B = abs(z2 - z1); 
            
            %   calculate gradient ::
            gradTh(iDepth) = (th2 - th1) / (z2 - z1); 
            
            %   calculate error ::
            deltaGradTh(iDepth) = sqrt(((deltaTh1 / B) ^ 2) ...
                                       + ((deltaTh2 / B) ^ 2)); 
            
        %   calculate for backward difference ::
        elseif iDepth == n
            
            %   get thorim values ::
            th1 = yth(iDepth - 1); 
            th2 = yth(iDepth); 
            
            %   get thorium error ::
            deltaTh1 = nth(iDepth - 1); 
            deltaTh2 = nth(iDepth); 
            
            %   get depth values ::
            z1 = x(iDepth - 1); 
            z2 = x(iDepth);
            
            %   calculate B ::
            B = abs(z2 - z1); 
            
            %   calculate gradient ::
            gradTh(iDepth) = (th2 - th1) / (z2 - z1); 
            
            %   calculate error ::
            deltaGradTh(iDepth) = sqrt(((deltaTh1 / B) ^ 2) ...
                                       + ((deltaTh2 / B) ^ 2)); 
            
        %   calculate for centered difference ::
        else
            
            %   get thorim values ::
            th1 = yth(iDepth - 1); 
            th2 = yth(iDepth + 1); 
            
            %   get thorium error ::
            deltaTh1 = nth(iDepth - 1); 
            deltaTh2 = nth(iDepth + 1); 
            
            %   get depth values ::
            z1 = x(iDepth - 1); 
            z2 = x(iDepth + 1);
            
            %   calculate B ::
            B = abs(z2 - z1); 
            
            %   calculate gradient ::
            gradTh(iDepth) = (th2 - th1) / (z2 - z1); 
            
            %   calculate error ::
            deltaGradTh(iDepth) = sqrt(((deltaTh1 / B) ^ 2) ...
                                       + ((deltaTh2 / B) ^ 2)); 
            
        end
        
    end
    %   spline smoothing approach removed: diagnostics confirmed smoothed ≈ raw FD within
    %   errorbars at all stations; residual term inflated σ at anchor points for no
    %   scientific gain on sparse profiles (8–15 obs per station) ::
    % deltaGradTh = deltaGradTh .^ 2;  % was variance for inverse-variance spline weights
    % weights  = 1 ./ deltaGradTh;
    % zBot     = min(rootsFt, ez);
    % W_anchor = max(weights);
    % x_aug    = [x;        mldVal;    zBot    ];
    % gTh_aug  = [gradTh';  0;         0       ];
    % wts_aug  = [weights'; W_anchor;  W_anchor];
    % [x_aug, sortIdx] = sort(x_aug);
    % gTh_aug  = gTh_aug(sortIdx);
    % wts_aug  = wts_aug(sortIdx);
    % gradThHat     = csaps(x_aug, gTh_aug, P_SPLINE_GRADIENT, x, wts_aug);
    % deltaGradThHat = sqrt((abs(gradThHat - gradTh') .^ 2) + deltaGradTh');

    %   hard-zero gradient and uncertainty outside active zone [mldVal, zBot] ::
    zBot = min(rootsFt, ez);
    gradTh(x < mldVal | x > zBot)      = 0;
    deltaGradTh(x < mldVal | x > zBot) = 0;

    %   store ::
    statIdx = find(gp15_obs.stationNo == sn & gp15_obs.depth <= maxDepth);
    gp15_grad.vertGrad(statIdx)      = gradTh;
    gp15_grad.vertGradError(statIdx) = deltaGradTh;
            
    %   clear arrays ::
    clear('gradTh', 'deltaGradTh'); 
            
end
close('all'); 

%%  end subroutine
