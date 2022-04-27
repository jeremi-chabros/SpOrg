clearvars; clc; close all;
load lookup_all.mat;
for ij = 1:max([S.const])
    % ij = 1;
    clearvars -except files ij

    construct_id = ij;
    q = [files.const];
    F = files(q == construct_id); % load just this construct

    % Here we obtain an appropriate/sorted file series
    T = struct2table(F);
    sortedT = sortrows(T, 'order');
    F = table2struct(sortedT);

    %% Set params
    fs = 25000;
    tolerance_ms = 5; % For artifacts
    bin_s = 20; % for moving average window
    method = 'thr3p5';
    spk_vec_all_recs = [];

    % Main loop
    for file = 1:length(F)

        % Load traces
        filename = F(file).name;
        load(filename)

        % Load spikes
        filename = convertStringsToChars(filename);
        filename = [filename(1:end-4) '.mat_spikes.mat'];
        load(filename);
        if file == 1
            threshold_list = thresholds;
        end

        % Filtering
        lowpass = 600;
        highpass = 8000;
        wn = [lowpass highpass] / (fs / 2);
        filterOrder = 3;
        [b, a] = butter(filterOrder, wn);
        filtered_data = filtfilt(b, a, double(dat));

        % Get some params & preallocate
        stim = files(file).stim;
        stim = [stim 15];
        duration_s = spikeDetectionResult.params.duration;
        duration_frames = ceil(fs*duration_s);
        rec_dur(file) = duration_s;

        num_chan = length(channels);
        spk_vec_all = zeros(1,  duration_frames);


        for j = 1:num_chan

            if ~ismember(channels(j), stim)
                spk_vec = zeros(1,  duration_frames);

                [stim_times, ~] = findStims(filtered_data(:,j));
                sps = spikeWaveforms{j}.(method);
                sps = sps(:,25);
                G = logical(and(sps<thresholds{j}.(method),sps>-50));

                spike_times = spikeTimes{j}.(method);
                spike_times = round(spike_times*25000);

                sp_times = spike_times(G);


                if ~isempty(stim_times)
                    artifact_locs = ismembertol(sp_times, stim_times,tolerance_ms*25,...
                        'DataScale', 1);
                    sp_times = sp_times(logical(~artifact_locs));
                end

                spk_times = sp_times;

                spk_vec(spk_times) = 1;
                spk_vec = spk_vec(1: duration_frames);

                %                 spikingFreq(j) = fs*sum(spk_vec)/duration_frames;
                spikingFreq(j) = length(spikeTimes{j}.(method))/duration_s;
                spk_vec_all = spk_vec_all+spk_vec;
            end
        end
        F(file).spikeFreq = spikingFreq;
        F(file).spikeFreqMean = sum(spikingFreq)/(60-(length(stim)));
        F(file).duration_frames = duration_frames;
        F(file).duration_s = duration_frames/fs;
        spk_vec_all_recs = horzcat(spk_vec_all_recs, spk_vec_all);

    end
    %% Get time points for labels
    for i = 1:length(rec_dur)
        rd = 0;
        for j = 1:i
            rd = rec_dur(1:j);
        end
        rdd(i) = sum(rd);
    end
    rdd = rdd-rdd(1);


    %% Plotting
    plot(movmean(spk_vec_all_recs, bin_s*fs),'k', 'linewidth', 2)

    hold on

    % Add lines with labels
    for i = 1:length(rdd)
        recname = F(i).name;
        recname = convertStringsToChars(recname);
        recname = recname(strfind(recname, 'Const')+18:end-4);
        xline(rdd(i)*fs,'r--',recname,'interpreter','none');
        hold on
    end

    % Aesthetics
    box off

    xlabel('Time (s)');
    xticks(linspace(1, rdd(end)*fs,10));
    xticklabels(round(linspace(1, rdd(end), 10)));

    ylabel('Relative spiking frequency');
    yticklabels(get(gca, 'ytick')*1000);

    set(gcf,'unit','normalized',...
        'innerposition',[1 1 1 1],...
        'color','w');

    title('210514_const1_cond5','interpreter','none'); % Could automate

    % Save

    plotPath = '/Users/jjc/mea/SpOrg/SpOrg_Plots';
    fname = F(1).name;
    fname = fname(1:end-4);
    fname = fname(1:strfind(fname, 'Condition')+9);

    if length(fname)>4
        exportgraphics(gcf, [plotPath filesep fname '.png'],'resolution',300);
        close(gcf)
        save(['/Users/jjc/mea/SpOrg/SpOrg_Results' filesep fname '_results.mat'], 'F', '-mat');
    else
        fname = F(1).name;
        fname = fname(1:end-4);
        surr = find(fname == '_', 1, 'last');
        fname = fname(1:surr-3);
        exportgraphics(gcf, [plotPath filesep fname '.png'],'resolution',300);
        close(gcf)
        save(['/Users/jjc/mea/SpOrg/SpOrg_Results' filesep fname '_results.mat'], 'F', '-mat');
    end
end
%%