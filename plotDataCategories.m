clearvars; clc;

% load('lookup_corrected1.mat');
load('lookup_corrected_thr3p5_2.mat');
S = files;
tags_ = {S.tag};


sgrp = [strcmp(tags_,'baseline')];
val_b = ([S(sgrp).spike_freq]);
val_b = rmoutliers([S(sgrp).spike_freq]);

sgrp = [strcmp(tags_,'stim')];
val_s = rmoutliers([S(sgrp).spike_freq]);
% val_s = val_s(val_s>3);

sgrp = [strcmp(tags_,'ttx')];
val_t = rmoutliers([S(sgrp).spike_freq]);

sgrp = [strcmp(tags_,'axotomy')];
val_a = rmoutliers([S(sgrp).spike_freq]);



notBoxPlot(val_b, 1);
hold on
notBoxPlot(val_s, 2);
hold on
notBoxPlot(val_t, 3);
hold on
notBoxPlot(val_a, 4);

h = gcf;
H = get(h);
set(gca, 'xtick', [1 2 3 4],...
    'xticklabels', {'baseline','stim','ttx','axotomy'});
xlabel('Recording condition')
ylabel('Spiking frequency [Hz]')
axis square


[p1,h,stats] = ranksum((val_b), (val_s));
[p2,h,stats] = ranksum((val_s), (val_t));
[p3,h,stats] = ranksum((val_t), (val_a));

hold on

sigstar({[1,2],[2,3],[3,4]},[p1, p2, p3])

exportgraphics(gcf, 'results.png', 'resolution', 300)


