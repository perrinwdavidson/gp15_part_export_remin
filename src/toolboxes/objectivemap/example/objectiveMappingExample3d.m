%%  objectiveMappingExample3d - objective mapping
%%-------------------------------------------------------------------------
%%  clear all
clear 
close all

%%  set data
%   data locations ::
Xvec = sort(10 * rand(10, 1)); 
Yvec = sort(10 * rand(10, 1));
Zvec = sort(10 * rand(10, 1));
[X, Y, Z] = ndgrid(Xvec, Yvec, Zvec); 

%   query grids ::
Xvecq = linspace(0, 10, 10)';
[Xq, Yq, Zq] = ndgrid(Xvecq, Xvecq, Xvecq); 

%   make samples ::
V = exp(-(((X - 5) .^ 2 + (Y - 5) .^ 2 + (Z - 5) .^ 2) ./ (2 * 1)));

%   flatten ::
Xin = [X(:), Y(:), Z(:)]; 
Xqin = [Xq(:), Yq(:), Zq(:)]; 
Vin = V(:); 
VVarin = zeros(size(Vin)); 

%   set max length scale ::
lenScaleBound = 1.2;

%%  objectively map data
[Vhat, VhatVar] = objectiveMapping(Xin, Vin, Xqin, lenScaleBound, VVarin);
Vhat = reshape(Vhat, 10, 10, 10); 
VhatVar = reshape(VhatVar, 10, 10, 10); 

%%  plot
figure; 
tl = tiledlayout(1, 3, 'tileSpacing', 'compact');
nexttile(); 
contourf(X(:, :, 5), Y(:, :, 5), V(:, :, 5)); 
colorbar(); 
title('Original');
nexttile(); 
contourf(Xq(:, :, 5), Yq(:, :, 5), Vhat(:, :, 5)); 
colorbar(); 
title('Objectively Mapped'); 
nexttile(); 
contourf(Xq(:, :, 5), Yq(:, :, 5), VhatVar(:, :, 5)); 
colorbar(); 
title('Mapping Standard Deviation'); 
set(gcf, 'position', [100 100 1300 400]); 

%%  end example
