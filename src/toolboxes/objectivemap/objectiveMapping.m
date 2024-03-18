function [VhatEstMean, VhatVar, lx] = objectiveMapping(X, V, Xq, lenScaleBound, VVar) 

    %%  normalize data
    %   calculate statistics ::
    Vmean = mean(V); 
    Vstd = std(V); 

    %   normalize ::
    Vnorm = (V - Vmean) ./ Vstd;

    %%  calculate initial coefficients
    %   set data variance ::
    varData = var(Vnorm); 

    %   make measurement variance matrix ::
    varMeasureData = diag(VVar); 

    %%  estimate length scale and variance 
    %   define function for correlation scale (theta) ::
    f = @(lenScale) -estLength(Vnorm, ...
                               (varData + varMeasureData) .* calcCovMat(squareform(pdist(X)), lenScale));  % 8.22 from Webster and Oliver (think nugget effect)

    %   find the length scale MLE:
    lenScaleMLE = fminbnd(f, 0, lenScaleBound);

    %   calculate the MLEs for sigma:
    [~, ~, varMLE] = estLength(Vnorm, ...
                               (varData + varMeasureData) .* calcCovMat(squareform(pdist(X)), lenScaleMLE));  

    %   set length scale ::
    lx = lenScaleMLE; %round(lenScaleMLE, 2);
    % disp(['Length scale = ' num2str(round(lx, 2))]);

    %   set variance ::
    varData = varMLE;
    C0 = varData;  % variance of the process
    
    %%  calculate noise variance (follow Bindoff and Wunsch 1992, Eq. (1))
    %   calculate square distance ::
    valSquare = pdist(Vnorm, 'euclidean') .^ 2;
    distSquare = pdist(X, 'euclidean');
    valSquare(distSquare > round(lx, 0)) = 0;
    
    %   calculate number of samples ::
    n = length(Vnorm); 
    
    %   calculate noise variance ::
    noiseVar = sum(valSquare, 'all') / (2 * n); 
    noiseVar = diag(noiseVar);

    %%  calculate data-data covariance matrix
    %  calculate distances ::
    rij = squareform(abs(pdist(X))); 
    expMat = exp(-(rij / lx) .^ 2); 

    %  calculate matrix ::
    Cdd = ((varData + varMeasureData) .* expMat) + noiseVar;

    %   calculate inverse ::
    CddInv = inv(Cdd);
    
    %%  calculate model-data covariance matrix
    %   set model variance ::
    varMod = varData;  % assume that variance is the same (think semivariogram, as opposed to Lerner et al., 2018)

    %   calculate distances ::
    rij = abs(pdist2(Xq, X)); 
    expMat = exp(-(rij / lx) .^ 2); 
    vij = abs(pdist2(zeros(size(Xq, 1), 1), VVar));

    %   calculate matrix ::
    Cmd = (varMod + vij) .* expMat;  % n.b.: we add in measurement error (think nugget)

    %%  calculate estimated mean
    estMean = sum(CddInv*Vnorm, 'all') / sum(CddInv, 'all');

    %%  calculate objective mapping (correcting following Deutsch, 1995)
    %%% calculate weights ::
    lmbda = Cmd * CddInv;
    lmbdaNew = lmbda;
    % lmbdaP = lmbda; 

    %%% calculate metrics ::
    % idxNeg = lmbda < 0;
    % lmbdaMean = mean(abs(lmbda(idxNeg))); 
    % cMean = mean(Cmd(idxNeg));

    %%% correct negative weights ::
    % lmbdaP(lmbda < 0) = 0;
    % lmbdaP((lmbda > 0) & (Cmd < cMean) & (lmbdaP < lmbdaMean)) = 0;
    % lmbdaNew = lmbdaP / sum(lmbdaP);

    %%% map ::
    VhatNorm = lmbdaNew * Vnorm;  
    VhatNormEstMean = lmbdaNew * (Vnorm - estMean) + estMean;

    %%  calculate sum of estimated mean (sanity check)
    sumEstMean = sum(CddInv*(Vnorm - estMean), 'all'); 
    
    %%  calculate model-model covariance matrix
    %   calculate model variance ::
    if size(Xq, 1) == 1

	% set model variance ::
	Cmm = C0;  % Eq. 8.35 from Webster and Oliver (*)

	% correct with measurement variance propagation ::
	Cmm = Cmm - ((lmbdaNew .^ 2) * VVar);  

    else

	% set model variance ::
    	varMod = var(VhatNormEstMean, 1);  % Bretherton, 1976 Eq. 10 OR Wunsch, 2006 Eq. 3.28 (**)

        % calculate distances ::
        rij = squareform(abs(pdist(Xq))); 
        expMat = exp(-(rij / lx) .^ 2); 

	% propagate measurement variance ::
	varProp = (lmbdaNew .^ 2) * VVar;

        % calculate matrix, correcting for measurement variance ::
        Cmm = (varMod - varProp) .* expMat; 

    end

    %%  calculate variance in estimate
    VhatVarNorm = Cmm - (lmbdaNew * Cmd');  % see (*) and (**) 

    %%  renormalize
    %   data ::
    Vhat = (VhatNorm .* Vstd) + Vmean;
    VhatEstMean = (VhatNormEstMean .* Vstd) + Vmean;
    
    %   variance ::
    %%% calculate ::
    VhatVar = diag((VhatVarNorm .* Vstd) + Vmean);

    %%% deal with negative variance ::
    VhatVar(VhatVar < 0) = 0;  % look up 'negative kriging weights' on how to handle differently, but this is a common way. due to close clustering of sample points around a query point.

end
