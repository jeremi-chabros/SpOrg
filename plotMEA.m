function [F, cbar] = plotMEA(X, A, Y, args)

% Description:
%   Plots custom heatmap with customizable elements:
%   - organization of the grid
%   - custom colormaps
%   - adjustable sizes and shapes of markers (faces + edges)
%   - option to input text into markers

% INPUT
%   X: vector of values for node colours
%   A: adjacency matrix specifying edge weights
%   Y: vector of values for node sizes
%   args: 'Name', Value pairs
%       'grd': vector of grounded (excluded) nodes; default: []
%       'c_map': (string) name of the custom heatmap, default: parula
%       'corners': (binary) flag to include corner coordinates; default: 0
%       'markerSize': size (area) of the markers; default: 2000
%       'marker': type of the markers; default 'o'
%       'edgeColor': color of the edges of valid nodes; default: 'none'
%       'cbarLimits':
%       'cbarTitle'
%       'cbar'
%       'edges'

% OUTPUT
%   F: handle to the figure
%   cbar: handle to the colorbar

% EXAMPLE
%   [F, cbar] = customHeatmap(X,...
%       'grd', [23 32],...
%       'corners', 1,...
%       'c_map', 'bamako',...
%       'markerSize', 2000,...
%       'marker', 'o',...
%       'edgeColor', 'none');
%
%    title('Spiking frequency heatmap',...
%        'fontsize', 14);
%    ylabel(cbar, {[''] ['Spiking frequency [Hz]']}, 'fontsize', 12);

% Author:
%   Jeremy Chabros, University of Cambridge, 2020
%   email: jjc80@cam.ac.uk
%   github.com/jeremi-chabros

%--------------------------------------------------------------------------
% CUSTOM COLORMAPS
%
% cmocean package
% cmocean colourmaps call master function with name as argument
% 'thermal','balance','haline','delta','solar','curl',
% 'ice','diff','gray','tarn','oxy','deep','dense','phase','algae',
% 'matter','turbid','topo','speed','amp','tempo','rain'
%
% See:
% https://uk.mathworks.com/matlabcentral/fileexchange/57773-cmocean-perceptually-uniform-colormaps
%
%--------------------------------------------------------------------------
% Scientific Colour Maps package
% Scientific colour maps are stored as .mat files
% addpath('/User/colourmaps/folder/')
%
% Available colourmaps:
% 'acton', 'bamako', 'batlow', 'berlin', 'bilbao', 'broc', 'brocO', 'buda',
% 'cork', 'corkO', 'davos', 'devon', 'grayC', 'hawaii', 'imola', 'lajolla',
% 'lapaz', 'lisbon', 'nuuk', 'oleron', 'oslo', 'roma', 'romaO', 'tofino',
% 'tokyo', 'turku', 'vik', 'vikO'
%
% See:
% https://zenodo.org/record/4153113#.X8j3gxP7TDI
% http://www.fabiocrameri.ch/colourmaps.php
%--------------------------------------------------------------------------

arguments
    X;
    A;
    Y;
    args.grd;
    args.c_map;
    args.corners;
    args.markerSize = 2000;
    args.marker = 'o';
    args.edgeColor = 'none';
    args.cbarLimits;
    args.cbarTitle;
    args.cbar;
    args.max;
    %     args.edges;
end

% Y = rescale(Y, 50,300);
A = rescale(A, 0, 4);

channels = [47;48;46;45;38;37;28;36;27;17;26;16;35;25;15;14;24;34;13;...
    23;12;22;33;21;32;31;44;43;41;42;52;51;53;54;61;62;71;63;...
    72;82;73;83;64;74;84;85;75;65;86;76;87;77;66;78;67;68;55;...
    56;58;57];

% Reference and grounded electrodes XY coords
ref = 15;

if isfield(args, 'grd')
    grd = args.grd;
    grd = intersect(grd,channels,'stable');
else
    grd = [];
end
% grd = channels(grd);

% if isfield(args, 'edges')
%     args.edges = [1,5];
% else

% Assign colours to values
% NOTE: ugly rescaling, probably needs a neater implementation
values = X;
values(61) = 0;
values(62) = args.max;
values = round(rescale(values, 1, 256));
values = values(1:60);

% List for the cmocean - calls master function
c_map_ocean = {'thermal','balance','haline','delta','solar','curl',...
    'ice','diff','gray','tarn','oxy','deep','dense','phase','algae',...,
    'matter','turbid','topo','speed','amp','tempo','rain', '-thermal',...
    '-balance','-haline','-delta','-solar','-curl','-ice','-diff','-gray',...
    '-tarn','-oxy','-deep','-dense','-phase','-algae','-matter','-turbid',...
    '-topo','-speed','-amp','-tempo','-rain'};

if isfield(args, 'c_map')
    
    c_map = args.c_map;
    
    if ismember(c_map, c_map_ocean)
        c_map = cmocean(c_map);
    else
        cmap_file = load(c_map);
        c_map = cmap_file.(c_map);
    end
    
    colormap(c_map);
else
    c_map = parula; % default colormap
end

if ~isfield(args, 'cbar')
    args.cbar=0;
    cbar=[];
end

%% Plot corner coordinates for reference
for i = 1:60
    nx(i) = channels(i);
    unitx(i) = rem(nx(i), 10);
    decx(i) = (nx(i)-unitx(i))/10;
end

% Plot edges
for i = 1:length(channels)
    for j=1:i
        if (channels(i)~=channels(j)) && (A(i,j)~=0)
            plot([decx(i), decx(j)], [unitx(i), unitx(j)],...
                'color', [.5 .5 .5 .5],...
                'linewidth', A(i,j));
        end
    end
    hold on;
end

% Plot corners if desired
if isfield(args, 'corners') && args.corners == 1
    textscatter(1, 8, string(18));
    hold on;
    textscatter(8, 1, string(81));
    hold on;
    textscatter(1, 1, string(11));
    hold on;
    textscatter(8, 8, string(88));
    hold on;
end


for i = 1:length(channels)
    
    n = channels(i);
    unit = rem(n, 10);
    dec = (n-unit)/10;
    
    % Specify what to do with ref electrode
    if n == ref
        
        f = scatter(dec, unit, args.markerSize, args.marker);
        % textscatter(dec, unit, "ref");
        f.MarkerFaceColor = 'none';
        f.MarkerEdgeColor = 'none';
        % f.MarkerEdgeColor = [0.5 0.5 0.5];
        
        % Specify what to do with grounded/stim electrode
    elseif ismember(n, grd)
        
        f = scatter(dec, unit, Y(i), args.marker);
        hold on;
        % textscatter(dec, unit, "ϟ", 'colordata', [1 1 1]);
        f.MarkerFaceColor = 'w';
        % f.MarkerEdgeColor = args.edgeColor;
        % f.MarkerFaceColor = c_map(values(i), :);
        f.LineWidth = 1;
        % f.MarkerFaceAlpha = 0.5;
        f.MarkerEdgeColor = 'k';
        
    else
        if Y(i) <= 50
            f = scatter(dec, unit, 30, 's');
        else
            f = scatter(dec, unit, Y(i), args.marker);
        end
        hold on;
        f.MarkerEdgeColor = args.edgeColor;
        f.LineWidth = 1;
        f.MarkerFaceColor = c_map(values(i), :);
    end
end
set(gca, 'ydir','reverse');
hold on;

% Set colorbar
colormap(c_map);

if isfield(args, 'cbarLimits') && numel(args.cbarLimits) == 2
    cbar_min = args.cbarLimits(1);
    cbar_max = args.cbarLimits(2);
else
    cbar_min = min(X);
    cbar_max = max(X);
end

if args.cbar
    
    cbar = colorbar;
    
    %     tks = linspace(cbar_min, cbar_max, 5);
    %     tks = round(tks, 3);
    %     tks1 = rescale(tks, 0, 1);
    %     cbar.Ticks = tks1;
    %     cbar.TickLabels = {round(tks(:),2)};
    cbar.Location = 'southoutside';
    cbar.TickDirection = 'out';
    cbar.Box = 'off';
    cbar.LineWidth = 1;
    if isfield(args, 'cbarTitle')
        set(get(cbar,'label'),'string',args.cbarTitle);
    end
    
    caxis([0 args.max]);
    %     cbar.Limits = [0 args.max];
    
    
else
    cbar = [];
end

% Set figure size and remove axes

xlim([0 9]);
ylim([0 9]);

set(gca,...
    'xcolor', 'none',...
    'ycolor', 'none',...
    'color', 'none');

% set(gcf, 'Position', [300 300 600 600]);
axis square
F = gcf;