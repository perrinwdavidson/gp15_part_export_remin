function [covMat] = calcCovMat(rij, ...
                               l)
%--------------------------------------------------------------------------
% About:
% Calculate the covariance function and corresponding matrix given the 
% inputted sample data. 
%
% Inputs:
% rij - distances between model and data points
% l - the length scale to be maximized.
%
% Outputs:
% covMat - correlation matrix
%--------------------------------------------------------------------------
    
    % find ratio of distances to correlation length scale:
    adjDist = rij / l; 

    % determine the appropriate correlation vector based on input: 
    covMat = exp(-(adjDist) .^ 2);
    
    % find size of correlation matrix:
    [m, ~] = size(covMat);
            
    % give same locations 1 covariance (diagonal):
    covMat(1 : m+1 : end) = 1;

end
