%%  calculate upwelling
%   find bounds of array ::
NUMLAT = length(u_mercator.latitude);
NUMLON = length(u_mercator.longitude);
NUMDEPTH = length(u_mercator.depth);
NUMTIME = length(u_mercator.time);

%   preallocate w array ::
w = NaN(size(u_mercator.uo));

%   get depth arrays ::
depth = u_mercator.depth;

%   iterate over all profiles ::
for iTime = 1 : 1 : NUMTIME
    for iLon = 1 : 1 : NUMLON
        for iLat = 1 : 1 : NUMLAT

            % get central velocity profiles ::
            uProfCentral = rmmissing(squeeze(u_mercator.uo(iLon, iLat, :, iTime)));
            vProfCentral = rmmissing(squeeze(u_mercator.uo(iLon, iLat, :, iTime)));
            uLonCentral = u_mercator.longitude(iLon);
            vLatCentral = v_mercator.latitude(iLat);
            
            % find length of profiles ::
            numSamp = min(length(uProfCentral), length(vProfCentral));
            
            % get forward velocity profiles and coordinates ::
            if iLon == NUMLON
                uProfForward = NaN(numSamp, 1);
                uLonForward = NaN;
            else
                uProfForward = rmmissing(squeeze(u_mercator.uo(iLon + 1, iLat, :, iTime)));
                uLonForward = u_mercator.longitude(iLon + 1);
            end
            if iLat == NUMLAT
                vProfForward = NaN(numSamp, 1);
                vLatForward = NaN;
            else
                vProfForward = rmmissing(squeeze(v_mercator.vo(iLon, iLat + 1, :, iTime)));
                vLatForward = v_mercator.latitude(iLat + 1);
            end    
            
            % get backward velocity profiles ::
            if iLon == 1
                uProfBackward = NaN(numSamp, 1);
                uLonBackward = NaN;
            else
                uProfBackward = rmmissing(squeeze(u_mercator.uo(iLon - 1, iLat, :, iTime)));
                uLonBackward = u_mercator.longitude(iLon - 1);
            end
            if iLat == 1
                vProfBackward = NaN(numSamp, 1);
                vLatBackward = NaN;
            else
                vProfBackward = rmmissing(squeeze(v_mercator.vo(iLon, iLat - 1, :, iTime)));
                vLatBackward = v_mercator.latitude(iLat - 1);
            end
            
            % skip if on land ::
            if isnan(numSamp)
                continue
            end
            
            % find length of all profiles ::
            numSamp = min([length(uProfCentral), length(vProfCentral), ...
                           length(uProfBackward), length(vProfBackward), ...
                           length(uProfForward), length(vProfForward)]);
            
            % specify arrays :: 
            uProfCentral = uProfCentral(1 : numSamp);
            vProfCentral = vProfCentral(1 : numSamp);
            uProfForward = uProfForward(1 : numSamp);
            vProfForward = vProfForward(1 : numSamp);
            uProfBackward = uProfBackward(1 : numSamp);
            vProfBackward = vProfBackward(1 : numSamp);
	    zProf = depth(1 : numSamp); 

            % calculate du/dx ::
            if iLon == 1  % forward finite difference
                
                % calculate du ::
                du = uProfForward - uProfCentral;  % m s-1
                
                % calculate dx ::
                dxDeg = uLonForward - uLonCentral;  % deg
                dx = repmat(deg2km(dxDeg, A), numSamp, 1);  % m
                
                % calculate du/dx ::
                dudx = du ./ dx; % sec-1
                
            elseif iLon == NUMLON  % backward finite difference
                
                % calculate du ::
                du = uProfCentral - uProfBackward;   % m sec-1
                
                % calculate dx ::
                dxDeg = uLonCentral - uLonBackward;  % deg
                dx = repmat(deg2km(dxDeg, A), numSamp, 1);  % m
                
                % calculate du/dx ::
                dudx = du ./ dx;  % sec-1
                
            else  % central difference
                
                % calculate du ::
                du = uProfForward - uProfBackward;  % m sec-1
                
                % calculate dx ::
                dxDeg = uLonForward - uLonBackward;  % deg
                dx = repmat(deg2km(dxDeg, A), numSamp, 1);  % m
                
                % calculate du/dx ::
                dudx = du ./ dx;  % sec-1
                
            end
            
            % calculate dv/dy
            if iLat == 1  % forward finite difference
                
                % calculate dv ::
                dv = vProfForward - vProfCentral;  % m sec-1
                
                % calculate dy ::
                dyDeg = vLatForward - vLatCentral;  % deg
                dy = repmat(deg2km(dyDeg, A), numSamp, 1);  % m
                
                % calculate du/dx ::
                dvdy = dv ./ dy;  % sec-1
                
            elseif iLat == NUMLAT  % backward finite difference
                
                % calculate dv ::
                dv = vProfCentral - vProfBackward;  % m sec-1
                
                % calculate dy ::
                dyDeg = vLatCentral - vLatBackward;  % deg
                dy = repmat(deg2km(dyDeg, A), numSamp, 1);  % m
                
                % calculate dv/dy ::
                dvdy = dv ./ dy;  % sec-1
                
            else  % central difference
                
                % calculate du ::
                dv = vProfForward - vProfBackward;  % m sec-1
                
                % calculate dx ::
                dyDeg = vLatForward - vLatBackward;  % deg
                dy = repmat(deg2km(dyDeg, A), numSamp, 1);  % m
                
                % calculate du/dx ::
                dvdy = dv ./ dy;  % sec-1
                
            end
            
            % calculate divergence per box ::
            horizDiv = - (dudx + dvdy);
            
            % integrate diverge for upwelling using trapezoidal integration ::
	    w(iLon, iLat, 1 : numSamp, iTime) = cumtrapz(zProf, horizDiv); 

        end
    end
    
    %   display time calculating ::
    whatTime = ['Finished calculating upwelling velocity for Time = ' num2str(iTime)];
    % disp(whatTime)
    
end

%%  by right hand rule, need to change sign (positive up) ::
w = -w; 

%%  end subroutine
