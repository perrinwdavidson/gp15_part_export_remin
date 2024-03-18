%%  objectiveMappingExample2d - objective mapping
%%-------------------------------------------------------------------------
%%  clear all
clear 
close all

%%  set data
%   data locations ::
Xvec = sort(10 * rand(10, 1)); 
Yvec = sort(10 * rand(10, 1));
[X, Y] = meshgrid(Xvec, Yvec); 

%   query grids ::
Xvecq = linspace(0, 10, 100)';
[Xq, Yq] = meshgrid(Xvecq, Xvecq); 

%   make samples ::
V = exp(-(((X - 5) .^ 2 + (Y - 5) .^ 2) ./ (2 * 1)));

%   flatten ::
Xin = [X(:), Y(:)]; 
Xqin = [Xq(:), Yq(:)]; 
Vin = V(:); 
VVarin = zeros(size(Vin)); 

%   set max length scale ::
lenScaleBound = 2;

%%  objectively map data
[Vhat, VhatVar] = objectiveMapping(Xin, Vin, Xqin, lenScaleBound, VVarin);

%%  plot
figure;
tl = tiledlayout(1, 3, 'tileSpacing', 'compact');
nexttile(); 
contourf(X, Y, V); 
colorbar(); 
title('Original');
nexttile(); 
contourf(Xq, Yq, reshape(Vhat, 100, 100)); 
colorbar(); 
title('Objectively Mapped'); 
nexttile(); 
contourf(Xq, Yq, reshape(sqrt(VhatVar), 100, 100)); 
colorbar(); 
title('Mapping Standard Deviation'); 
set(gcf, 'position', [100 100 1300 400]); 

%%  end example
