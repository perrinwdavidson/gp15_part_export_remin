function [sSquared] = calcVar(vi, ...
                              di, ...
                              lenScale)
  
    %   calculate distances ::
    modeDists = pdist(di);
    
    %   calculate differences ::
    modeDiffs = (pdist(vi) .^ 2) / (2 * length(vi));
    
    %   specify distances ::
    modeDiffs(modeDists > lenScale) = NaN;
    modeDiffs = rmmissing(modeDiffs); 
    
    %   calculate variance ::
    sSquared = sum(modeDiffs, 'all'); 
                          
end

