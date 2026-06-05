%%  objectiveMappingExample - objective mapping
%%-------------------------------------------------------------------------
%%  clear all
clear 
close all
warning('off');

%%  set data
%   data locations ::
X = [1, 2, 3, 4, 5, 6, 7]';
V = [5, 6, 7, 6, 4, 5, 4]';
VVar = zeros(size(V)); 

%   query grids ::
Xq = linspace(0, max(X), 1000)';

%   set max length scale ::
lenScaleBound = 50;

%%  objectively map data
[Vhat, VhatVar] = objectiveMapping(X, V, Xq, lenScaleBound, VVar);

%%  plot
figure;
    %   set arrays ::
    Xplot = [Xq; flipud(Xq)];
    Yplot = [Vhat - sqrt(VhatVar); flipud(Vhat + sqrt(VhatVar))];
    hold on
    scatter(X, V, 'markerEdgeColor', [0 0 1], 'lineWidth', 1.5, ...
            'marker', 'o', 'sizeData', 50, 'lineWidth', 1.5, ...
            'markerFaceColor', [1, 1, 1]);
    plot(Xq, Vhat, 'color', [1 0 0], 'lineWidth', 2)
    h = fill(Xplot', Yplot', [0.8, 0.8, 0.8]);
    set(h, 'edgeColor', 'none')
    hold off
    grid on
    box on
        title('Objectively Mapped Data with Standard Deviation')
        legend('Sample', 'Mapped', 'Standard Deviation')
            chi = get(gca, 'children');
            set(gca, 'children',flipud(chi))

%%  end example
