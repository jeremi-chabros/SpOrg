clearvars; clc;
lkp_name = 'lookup_all';
load(lkp_name);

files = S;

method = 'thr3p5';
min_freq_Hz = 1;
filelist = {'190830_slice1','190830_slice2', '210514_Const1', '201029_Const3', '201029_Const2', '201029_Const1'};
suffix = '_sporg_e.mat';
% lower_threshold = -20;

% Need to get first baseline recs
T = struct2table(files); % convert the struct array to a table
sortedT = sortrows(T, 'order'); % sort the table by 'order'
files = table2struct(sortedT); % change it back to struct array if necessary

chan_list = [47;48;46;45;38;37;28;36;27;17;26;16;35;25;15;14;24;34;13;23;12;22;33;21;32;31;44;43;41;42;52;51;53;54;61;62;71;63;72;82;73;83;64;74;84;85;75;65;86;76;87;77;66;78;67;68;55;56;58;57];
%%
for n = 1:length(files)
    % for n = 7:10
    clear lower_threshold stim_interval;
    filename = files(n).name;
    if startsWith(filename, filelist)

        load(filename);
        load([filename(1:end-4) '.mat_spikes.mat']);
        ordinal = find(strcmp(filelist,filename(1:13)));
        load([filelist{ordinal} suffix]);

%         if ismember(n, 8:14)
%             thr_id{n-7} = spikeDetectionResult.params.mad; % the value itself is +ve
%         elseif n <= 7
%             thr_id{n} = spikeDetectionResult.params.mad;
%         end

        if ismember(n, 1:7)
            thr_id{n} = spikeDetectionResult.params.mad;
        end

        % thr_id is a cell {no_constructs, 1} containing [1,60] vector of
        % thresholds
        % so for a given recording thr_id{file(n).const} contains all
        % thresholds and thr_id{file(n).const(i) contains the value for
        % particular electrode

        stim = files(n).stim;
        stim = [stim 15];
        [~,locs] = ismember(stim, channels);
        locs = [locs, 15];

        lowpass = 600;
        highpass = 8000;
        wn = [lowpass highpass] / (fs / 2);
        filterOrder = 3;
        [b, a] = butter(filterOrder, wn);
        filtered_data = filtfilt(b, a, double(dat));

        for i = 1:length(channels)

            spike_freq(i) = 0;
            chan = channels(i);
            ypos = rem(chan,10);
            xpos = (chan-ypos)/10;

            indx = sub2ind([8,8],xpos,ypos);

            if ~ismember(chan, stim) && ismember(chan, electrode_id)

                [stim_times, stim_interval] = findStims(filtered_data(:,i));
                sps = spikeWaveforms{i}.(method);
                sps = sps(:,25);

                %             if strcmp(S(n).tag, 'baseline') || strcmp(S(n).tag, 'stim')
                %                 lower_threshold = -spikeDetectionResult.params.mad(i);
                %             else
                %                 lower_threshold = -20;
                %             end

                lower_threshold = -thr_id{files(n).const}(i);

                G = logical(and(sps<lower_threshold,sps>-100));

                spike_times = spikeTimes{i}.(method);
                spike_times = round(spike_times*25000);

                sp_times = spike_times(G);
                tolerance_ms = 5;

                if ~isempty(stim_times)
                    artifact_locs = ismembertol(sp_times, stim_times,tolerance_ms*25,...
                        'DataScale', 1);
                    sp_times = sp_times(logical(~artifact_locs));
                end

                spike_freq(i) = length(sp_times)/spikeDetectionResult.params.duration;

            end
            clear lower_threshold;
        end

        num_active = sum(spike_freq>min_freq_Hz);
        num_active = length(electrode_id);
        files(n).no_active_el = num_active;
        files(n).spike_freq = sum(spike_freq)/(60-length(locs));
        files(n).spike_freq = spike_freq;
%         files(n).spike_freq_alt = sum(spike_freq)/num_active;
        files(n).spike_freq_alt = sum(spike_freq)/length(electrode_id);
        files(n).duration_s = spikeDetectionResult.params.duration;


        if stim_interval
            files(n).stim_freq = fs*1./(mean(stim_interval));
        else
            files(n).stim_freq = 0;
        end


        files(1).volts = thr_id;
        save([lkp_name(1:end-4) 'thr3p5_sporg_only.mat'],'files');
    end
end