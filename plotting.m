clearvars; clc;
lkp_name = 'lookup_final.mat';
load(lkp_name);
spkmet = 'mea';
chan_list = [47;48;46;45;38;37;28;36;27;17;26;16;35;25;15;14;24;34;13;23;12;22;33;21;32;31;44;43;41;42;52;51;53;54;61;62;71;63;72;82;73;83;64;74;84;85;75;65;86;76;87;77;66;78;67;68;55;56;58;57];
%%
for n = 1:length(files)
% for n = 28:28
    close all
    filename = files(n).name;
    load(filename);

    load([filename(1:end-4) '.mat_spikes.mat']);
    
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
    
    % Plot grid
    clc
    stim = files(n).stim;
    sporg = files(n).sporg;
    figure
    
    for i = 1:length(channels)
        spike_freq(i) = 0;
        chan = channels(i);
        ypos = rem(chan,10);
        xpos = (chan-ypos)/10;
        
        indx = sub2ind([8,8],xpos,ypos);
        
        subplot(8,8,indx)
        plot(filtered_data(:,i), 'color','k')
        
        if ismember(chan, sporg)
            
            [stim_times, ~] = findStims(filtered_data(:,i));
            sps = spikeWaveforms{i}.(spkmet);
            sps = sps(:,25);
            G = logical(and(sps <- 0,sps>-100));
            
            spike_times = spikeTimes{i}.(spkmet);
            spike_times = round(spike_times*25000);
            
            sp_times = spike_times(G);
            tolerance_ms = 5;
            
            if ~isempty(stim_times)
                artifact_locs = ismembertol(sp_times, stim_times,tolerance_ms*25,...
                    'DataScale', 1);
                sp_times = sp_times(logical(~artifact_locs));
            end
            
            spike_freq(i) = length(sp_times)/spikeDetectionResult.params.duration;
            
            
            hold on
            scatter(sp_times, repmat(25,[length(sp_times),1]), 20,...
                'v', 'filled',...
                'markerfacecolor', [255, 69, 0]/256,...
                'markeredgecolor', 'k');
            
            if length(stim_times) > 1
                hold on
                scatter(stim_times, repmat(50,[length(stim_times),1]), 20,...
                    's', 'filled',...
                    'markerfacecolor', [5, 92, 157]/256,...
                    'markeredgecolor', 'none');
            end
        end
        
        xlim([1 25000*12])
        ylim([-50 50])
        
        
        if ismember(chan, stim)
            set(gca,'xcolor','r','ycolor','r',...
                'xtick',[],'ytick',[],...
                'linewidth', 1);
        else
            box off
            axis off
        end
        axx(i)=gca;
        
    end
    subplot(8,8,64)
    plot(filtered_data(1,i), 'color','w')
    
    xlim([1 25000*12])
    xlabel('Time (s)')
    xt = get(gca, 'XTick');
    
    ylim([-50 50])
    ylabel('Voltage (uV)')
    
    set(gca, 'XTick', [0,xt(end)], 'XTickLabel', [0,xt(end)]/25000,'tickdir','out');
    
    box off
    files(n).spike_freq = sum(spike_freq)/(60-length(locs));
    % linkaxes(axx,'y')
    sgtitle(files(n).name, 'interpreter','none','fontsize', 24)
    set(gcf,'color','w','unit','normalized','outerposition',[0 0 .625 1]);
    exportgraphics(gcf, ['/Users/jjc/mea/SpOrg/SpOrg_Plots/' filename(1:end-4) '.png'], 'resolution', 300);
    
    % Plot heatmap
    close all
    [F, cbar] = plotMEA(spike_freq, zeros(60), 800*ones(60,1),...
        'grd', files(n).stim,...
        'cbar', 1,...
        'c_map', 'thermal',...
        'corners',0,...
        'max', 20,...
        'cbarTitle', 'Spiking frequency (Hz)');
    set(gcf,'color','w');
    
    title(filename(1:end-4),'interpreter','none');
    
    save([lkp_name(1:end-4) '_corrected.mat'],'files');
    exportgraphics(gcf, ['/Users/jjc/mea/SpOrg/SpOrg_Plots/' filename '.png'],'resolution',300);
end


%%
figure
trace = filtered_data(:,6);
plot(trace,'k','linewidth',1)
% xlim([25000*10 25000*10.05])
ylim([-50 50])

xlabel('Time (ms)','fontsize', 14)
ylabel('Voltage (uV)','fontsize', 14)
box off

xt = get(gca, 'XTick');
set(gca, 'XTick', xt, 'XTickLabel', xt/25,'tickdir','out');

% title(['channel ' num2str(channels(i))])

set(gcf,'color','w','unit','normalized','outerposition',[0 0 1 1]);


%%
j = 21;
close all
tl = tiledlayout(5,1,'tilespacing', 'none', 'padding','none');
title(tl, files(n).name, 'interpreter', 'none','fontsize', 14)
xlabel(tl,'Time (ms)','fontsize', 14)
ylabel(tl,'Voltage (uV)','fontsize', 14)

for i = j:(j+4)
    nexttile
    %     plot(dat(:,i)-mean(dat(:,i)),'k','linewidth',1)
    plot(filtered_data(:,i),'k','linewidth',1)
    xlim([25000 25000*2])
    %     ylim([-100 100])
    box off
    
    xt = get(gca, 'XTick');
    set(gca, 'XTick', xt, 'XTickLabel', xt/25,'tickdir','out')
    title(['channel ' num2str(channels(i))])
    
    if i~=j+4
        set(gca,'xcolor','none');
    end
end

linkaxes(tl.Children(:),'y')
set(gcf,'color','w','unit','normalized','outerposition',[0 0 1 1]);



