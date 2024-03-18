%% olsDataExample -- performing ols on linear data 
%% --------------------------------------------------------
%  configure environment ::
close('all'); 
clear; 

%  make data ::
x = []; 
y = []; 

% regress dasta ::
[beta0Hat, beta1Hat, q, ci] = ols(x, y, 1, 1);
