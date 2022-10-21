clearvars; clc;

% load('lookup_corrected1.mat');
load('lookup_corrected_thr3p0_2.mat');
S = files;
tags_ = {S.tag};

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

for i = 1:max([S.const])

    surr_var = ([S.const]==1);
    surr = max([S(surr_var).order]);


    for j = surr:-1:1
        lh = plot(val.(tagnames{1}), 's-', 'color', [.5 .5 .5]);
        lh.Color = [lh.Color .5];
        hold on;
    end
end
