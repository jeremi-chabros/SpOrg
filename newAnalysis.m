clearvars; clc;
load("lookup_final.mat");
f = readtable('SpOrgs.csv');
method = 'thr3p0';
multiplier = 3.5;

% There are multiplicates of stim cond for const2
files(6).order = 2;
files(7).order = 2;
%%
% NOTES:
% - obtain MAD from baseline rec and use thr based on that in subsequent
% recordings
% - ONLY use the electrodes at SpOrg side
% - TODO: tragically not optimised with no pre-allocations and loads of
% overwriting variables

for i = 1:length(files)

    filename = files(i).name;
    sporg_el = files(i).sporg;

    load(filename); % Load raw data
    % Filtering raw data for artifact removal
    lowpass = 600;
    highpass = 8000;
    wn = [lowpass highpass] / (fs / 2);
    filterOrder = 3;
    [b, a] = butter(filterOrder, wn);
    filtered_data = filtfilt(b, a, double(dat));

    % Load spike detection results
    load([filename, '_spikes.mat']);

    sporg_el(sporg_el==85) = []; % remove the noisy electrode 85

    enum = 0;
    spk_fr_by_el = zeros(1,60);
    clear spikefreq
    for j = sporg_el % Run only for SpOrg electrodes
        enum = enum + 1;
        elcoord = find(channels == j);

        if files(i).order == 1
            thr = thresholds{elcoord}.('thr3p0')/3.0; % Extract 'raw' thr from multiplier
        end

        [stim_times, stim_interval] = findStims(filtered_data(:,elcoord));
        sps = spikeWaveforms{elcoord}.(method);
        sps = sps(:,25); % waveform peak
        G = logical(and(sps < multiplier*thr, sps > -70));

        sp_times = spikeTimes{elcoord}.(method);
        spike_times = round(sp_times*fs);

        sp_times = spike_times(G);
        tolerance_ms = 5;

        if ~isempty(stim_times)
            artifact_locs = ismembertol(sp_times, stim_times,tolerance_ms*fs/1000,...
                'DataScale', 1);
            sp_times = sp_times(logical(~artifact_locs));
        end
        spk_filtered = sp_times;
        numspikes(enum) = length(spk_filtered);
        spikefreq(enum) = numspikes(enum)/spikeDetectionResult.params.duration;
        spk_fr_by_el(elcoord) = numspikes(enum)/spikeDetectionResult.params.duration;

    end
    files(i).spikefreq = median(spikefreq); % get median
    files(i).spikeFrCh = spikefreq; % get spiking frequency by channel
end

%%
for i = 1:max([files.const])
    constf = files([files.const] == i);
    t = struct2table(constf);
    s = sortrows(t, 'order');
    constf = table2struct(s);
    cntr = 2;
    surr(i,1) = 1;
    for j = 1:max([constf.order])-1
        surr(i,cntr) = 1+((constf(j+1).spikefreq - constf(1).spikefreq)/constf(1).spikefreq);
        cntr = cntr+1;
    end
end

save('results.mat', 'surr', 'files');
%%
load('results.mat')

SEM = nanstd(surr)/sqrt(length(surr));

err = repmat(SEM,[5,1]);
bar(nanmedian(surr-1))

hold on
x = 1:5;
data = nanmedian((surr-1))';
errhigh = SEM;
errlow  = SEM;

% bar(x,data, 'facealpha', .5, 'edgecolor', 'none')             

% hold on

er = errorbar(x,data,errlow,errhigh);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  

hold on

surr(surr==0)=nan;
scatter(1:5,surr-1, 50, 'fill', 'markeredgecolor', 'k')

xlabel('Condition')
xticks = 1:4;
xticklabels({'Baseline','Stim','Post-stim','TTX', 'Axotomy'});
ylabel('Firing rate (% baseline)')
yticklabels(yticks*100)
set(gcf, 'color', 'w')
set(gca, 'tickdir', 'out')
box off
axis padded
legend({'','',files([files.order] == 1).name}, 'interpreter', 'none')

%%
surr(surr==0) = nan;
scatter(2:5,surr(:,2:end),'fill')
hold on
plot(surr','--')
hold on
% boxplot(surr)
axis padded
box off
set(gcf, 'color', 'white')
xlabel('Condition')
% xticks = 2:5:3;
xticklabels({'Baseline', '', 'Stim','', 'Post-stim','', 'TTX', '', 'Axotomy'});
% xticklabels({'Baseline','Stim','Post-stim','TTX', 'Axotomy'});
ylabel('Firing rate (% baseline)')
yticklabels(yticks*100)
set(gcf, 'color', 'w')
set(gca, 'tickdir', 'out')
legend({files([files.order] == 1).name}, 'interpreter', 'none')

%%

% values needs to be vector of size [60x1]
for i = 1:length(files)
% for i = 15:15
    filename = files(i).name;

    
    channels = [47;48;46;45;38;37;28;36;27;17;26;16;35;25;15;14;24;34;13;...
    23;12;22;33;21;32;31;44;43;41;42;52;51;53;54;61;62;71;63;...
    72;82;73;83;64;74;84;85;75;65;86;76;87;77;66;78;67;68;55;...
    56;58;57];

    spike_freq = zeros(1,60);
    files(i).sporg(files(i).sporg==85) = [];
    sporg_els = ismember(channels, files(i).sporg);
    sporg_els = sporg_els(sporg_els ~= 85);
    spike_freq(sporg_els) = files(i).spikeFrCh;

    grd = ~ismember(channels, files(i).sporg);
    grd = channels(grd == 1);
    grd = vertcat(grd, (files(i).stim)');
    grd = [grd; 85];

    [F, cbar] = plotMEA(spike_freq, zeros(60), 800*ones(60,1),...
        'grd', grd,...
        'cbar', 1,...
        'c_map', 'thermal',...
        'corners',0,...
        'max', 20,...
        'cbarTitle', 'Spiking frequency (Hz)');
    set(gcf,'color','w');

    title(filename(1:end-4),'interpreter','none');
    exportgraphics(gcf, ['/Users/jjc/mea/SpOrg/Plots/' filename(1:end-4) '.png'], 'resolution', 300);
    close all
end

%%
allvec = [];
Ords = [];
for o = 1:4
    hord = [files.order];
    hord = hord==o;
    newvec = files(hord).spikeFrCh;
    allvec = [allvec newvec];
    ords = zeros(length(newvec),1)+o;
    Ords = [Ords; ords]
end

violinplot(allvec,Ords);
set(gcf,'color','w')

box off
axis padded
xlabel("condition")
xticks([1,2,3,4])
xticklabels(["baseline","stim","post-stim","axotomy"]);
ylabel("Spiking frequency [Hz]")
set(gcf,'color','w')
