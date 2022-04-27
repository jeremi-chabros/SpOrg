clearvars; clc;
load('lookup.mat');

f = dir('/Users/jjc/mea/SpOrg/SpOrg_Spikes/*20*.mat');
for n = 1:length(f)
    close all
    load(f(n).name);

    %%
    method = 'thr3p5';
    stim = files(n).stim;
    [~,locs] = ismember(stim, channels);
    locs = [locs, 15];

    for i = 1:length(spikeWaveforms)

        sps = spikeWaveforms{i}.('thr3p5');
        sps = sps(:,25);
        if ~ismember(i,locs)
            count_abs(i) = sum(and(sps<-25,sps>-50))/spikeDetectionResult.params.duration;
            count_thr(i) = length(spikeTimes{i}.(method))/spikeDetectionResult.params.duration;
        else
            count_abs(i)=0;
        end
    end
    files(n).count_abs = sum(count_abs)/(60-length(locs));
    files(n).count_thr = sum(count_thr)/(60-length(locs));
    %%
    % close all

    [F, cbar] = plotMEA(count_abs, zeros(60), 800*ones(60,1),...
        'grd', files(n).stim,...
        'cbar', 1,...
        'c_map', 'thermal',...
        'corners',0,...
        'max', 20);
    set(gcf,'color','w');


plotPath = '/Users/jjc/mea/SpOrg/SpOrg_Plots';
    fname = f(n).name;
    fname = fname(1:strfind(fname, '.mat_spikes.mat')-1);

    title(fname,'interpreter','none');

    % save('../lookup.mat','files');
    exportgraphics(gcf, [plotPath filesep fname '.png'],'resolution',300);
end