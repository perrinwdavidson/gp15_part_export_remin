function [lVal, muPartial, varPartial] = estLength(Vnorm, ...
                                                   covMat)                                            
%--------------------------------------------------------------------------
% About:
% Computes the likelihood of multivariate normal random variables 
% for a given correlation matrix.
%
% Input:
% Vnorm - nx by nSamples matrix of realisations of a Gaussian process
% covMat - correlation matrix
%
% Output:
% lVal - likelihood value from evaluation of the likelihood function
% muPartial - MLE value of mu for given correlation matrix
% varPartial - MLE value of variance for given correlation matrix
%--------------------------------------------------------------------------

    % find dimensions of values at sample points:
    [nx, nSamples] = size(Vnorm);

    % get corrMat solution:
    oneVec = ones(nx, 1);
    oneVecLarge = repmat(oneVec, nSamples, 1);
    
    % only use 1d array inputs:
    if nSamples > 1
        disp('Needs to be 1d array.')
        return
    end

    % solve equation \rho*s = 1 from (33b) of [1]:
    s = covMat \ oneVec;
    
    % put s into matrix form:
    sLarge = repmat(s, nSamples, 1);
    
    % preallocation array for solution to \rho*r = x from (33a) of [1]:
    rLarge = zeros(nx * nSamples, 1);

    % for each vector within the samples:
    for iSample = 1 : nSamples 
        
        % solve (33a) per vector:
        r = covMat \ Vnorm(:, iSample);
        
        % get appropriate indices and put solution in r matrix:
        ind2 = nx * iSample;
        ind1 = ind2 - (nx - 1);
        rLarge(ind1:ind2) = r;
        
    end
    
    % put samples into vector:
    flatSamples = Vnorm; %reshape(Vnorm, [numel(Vnorm), 1]);

    % calculate determinant of correlation matrix ^[1]:
    logDet = (nSamples * 2) * sum(log(diag(chol(covMat))));
    
    % calculate the mean MLE, mu:
    muPartial = (oneVecLarge' * rLarge) / (oneVecLarge' * sLarge);
    
    % calculate the variance MLE, sigma^2:
    varPartial = (1 / (nx*nSamples)) ...
               * (flatSamples - (muPartial * oneVecLarge))' ...
               * rLarge;

    % calculate likelihood of the interpolation sample given \phi = {muMLE,
    % varMLE, \theta}, for which \theta (the correlation length) is 
    % presribed here:
    lVal = (-((nx*nSamples) / 2)) * log(varPartial) - (0.5*logDet);
    
    % calculate standard deviation:
    sigmaPartial = sqrt(varPartial);

end