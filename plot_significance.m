clearvars; clc;

% load('lookup_corrected1.mat');
load('lookupthr3p0_sporg_only.mat');
S = files;
tags_ = {S.tag};

% tagnames = unique(tags_);
tagnames = {'baseline', 'stim', 'ttx', 'axotomy'};
rmout = 0;

if rmout
    for i = 1:length(tagnames)
        sgrp = [strcmp(tags_, tagnames{i})];
        val.(tagnames{i}) = rmoutliers([S(sgrp).spike_freq_alt]);
        hold on
    end
else
    for i = 1:length(tagnames)
        sgrp = [strcmp(tags_, tagnames{i})];
        val.(tagnames{i}) = ([S(sgrp).spike_freq_alt]);
        hold on
    end
end


%
clc;
bar_width = 0.3;
jitWidth = bar_width;
sep = 3/2*bar_width;

pos = 1:4;
pos = 3/2*pos;
pos = [pos-sep; pos; pos+sep];

% colors = linspecer(3);
colors = [[0, 113, 188]/256; [182, 50, 28]/256; [0,0,0]; [30,76,16]/256];
mrks = {'o', '^', 's', 'd'};

for i = 1:numel(tagnames)
    pp = pos';
    boxplot(val.(tagnames{i})',...
        'widths', bar_width,...
        'colors', colors(i,:),...
        'symbol','',...
        'positions', pp(i));
    set(findobj(gca,'type','line'),'linew',1)
    hold on
    for j = 1:length(pos)
        jitter = linspace(0, bar_width, length(val.(tagnames{i})));
        s = scatter(pp(i)-jitWidth/2+jitter, (val.(tagnames{i})));
        s.Marker = mrks{i};
        s.MarkerFaceColor = colors(i,:);
        s.MarkerEdgeColor = colors(i,:);
        s.MarkerEdgeColor = 'w';
        s.SizeData = 50;
        s.MarkerFaceAlpha = 0.3;
    end
end

h = gcf;
H = get(h);
pos_2 = 1:4;

set(gca, 'xtick', pos(1,:),...
    'xticklabels', tagnames);

xlabel('Recording condition')
ylabel('Spiking frequency [Hz]')
axis square

val_perms = nchoosek(1:4,2);

for i = 1:length(val_perms)
    valname_1 = tagnames{val_perms(i,1)};
    valname_2 = tagnames{val_perms(i,2)};
    [p(i), h, stats] = ranksum(val.(valname_1), val.(valname_2));
end
% Phantom plots for the legend
h1=plot(nan, 'o-', 'color', colors(1,:), 'linewidth', 1,...
    'markerFaceColor', colors(1,:));
hold on
h2=plot(nan, '^-',  'color', colors(2,:), 'linewidth', 1,...
    'markerfacecolor', colors(2,:));
hold on
h3=plot(nan, 's-',  'color', colors(3,:), 'linewidth', 1,...
    'markerfacecolor', colors(3,:));
h4=plot(nan, 'd-',  'color', colors(4,:), 'linewidth', 1,...
    'markerfacecolor', colors(4,:));


hold on
axis padded
argpos = pos(2,:);
argpos = argpos(val_perms)-(1.5*bar_width);
sigstar({argpos(1,:),argpos(2,:),argpos(3,:),argpos(4,:),argpos(5,:),argpos(6,:)},[p])
box off
set(gcf, 'color','w', 'position', [300,300,1000,1000])

% Legend
l = legend([h1 h2 h3 h4]);
l.Title.String = 'Condition:';
l.String = tagnames;
l.Location = 'northeastoutside';
l.Box = 'off';
ax=gca;
ax.LineWidth = 1;

exportgraphics(gcf, 'results.png', 'resolution', 300)