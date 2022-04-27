clearvars; clc;
recs = dir('*-0.3*.mat');
recs = rmfield(recs,'folder');
recs = rmfield(recs,'date');
recs = rmfield(recs,'bytes');
recs = rmfield(recs,'isdir');
recs = rmfield(recs,'datenum');
[~,index] = sortrows({recs.name}.'); recs = recs(index); clear index
%%
fr = struct('Baseline',[],'Stim',[],'PostStim',[]);
max_count = 0;
method = 'thr4p5';
for q = 1:length(recs)
% for q = 1:8
    recname = recs(q).name;
    load(recname);
    for i = 1:60
        try
            counts(i) = length(spikeTimes{i}.(method))/spikeDetectionResult.params.duration;
        catch
            counts(i) = 0;
        end
    end
    
    if max_count <= max(counts)
        max_count = max(counts);
    end
    
    vals_sorted = sort(counts, 'descend');
    mean_fr(q) = mean(vals_sorted(1:50));
%     n_active = sum(counts>1);
    
    
    if contains(recname, 'Baseline')
        recs(q).cond = 'Baseline';
        fr.Baseline(end+1) = mean_fr(q);
    elseif contains(recname, 'PostStim')
        recs(q).cond = 'PostStim';
        fr.PostStim(end+1) = mean_fr(q);
    else
        recs(q).cond = 'Stim';
        fr.Stim(end+1) = mean_fr(q);
    end
    
    recs(q).spfr = mean_fr(q);
    recname = strrep(recs(q).name,'_', ' ');
    recname = strrep(recname,'p', '.');
    recname = recname(1:strfind(recname, 'L')-2);
    recs(q).name = recname;
    
    [~,cbar] = plotMEA(counts, zeros(60,60), ones(1,60)*750,...
    'c_map', 'thermal',...
    'cbar', 1,...
    'cbarTitle', 'Spiking frequency (Hz)',...
    'corners', 0,...
    'grd', [],...
    'cbarLimits',[0 max_count]);

axis square

set(gcf, 'color', 'w');
set(gcf,'units','pixels');
set(gca,'units','pixels');
w_pos = get(gcf, 'position');
% set(gca, 'position', [0 0 w_pos(3) w_pos(4)]);

% exportgraphics(gcf, ['/Users/jeremi/Desktop/Hmaps/' recname '.png'], 'resolution',300);
%     recs(q).n_active = n_active;
close all;
    
end
%%
n = 10;
load(recs(n).name);
close all;

% stim = [16 17 26 27 28 36 37 38 46 47 48];
r = [];
% for j = 1:length(stim)
%   r(j) = find(channels == stim(j));
% %   color_data(r(j),:) = [134 156 211]/255;
% end

for i = 1:length(spikeTimes)
    try
        counts(i) = length(spikeTimes{i}.(method))/spikeDetectionResult.params.duration;
    catch
        counts(i) = 0;
    end
end

[~,cbar] = plotMEA(counts, zeros(60,60), ones(1,60)*750,...
    'c_map', 'thermal',...
    'cbar', 0,...
    'cbarTitle', 'Spiking frequency (Hz)',...
    'corners', 0,...
    'grd', r,...
    'cbarLimits',[0 max_count]);

cb = colorbar('location','southoutside');
cb.Box = 'off';
cb.TickDirection = 'out';
cb.LineWidth = 1;
% cb.XLim = [0 max_count];
caxis([0 max(counts)])

axis square

set(gcf, 'color', 'w');
set(gcf,'units','pixels');
set(gca,'units','pixels');
w_pos = get(gcf, 'position');
% set(gca, 'position', [0 0 w_pos(3) w_pos(4)]);

recname = strrep(recs(n).name,'_', ' ');
recname = strrep(recname,'p', '.');
recname = recname(1:strfind(recname, 'L')-2);
title(recname)


%%
 [~,index] = sortrows({recs.cond}.'); recs = recs(index); clear index;
 a = [recs(1:13).spfr];
 qqplot(a);


%%
a = [];
% a = zeros(1,9);
a(1:4) = [recs(13:16).spfr];
a(5:8) = [recs(17:20).spfr];

close all;
clr = zeros(3,8);
clr(:,1:4) = repmat([223 145 167]/255,4,1)';
clr(:,5:8) = repmat([134 156 211]/255,4,1)';
% clr(:,7:end) = repmat([178 198 65]/255,2,1)';
b = bar(a,'facecolor', 'flat');
b.CData = clr';
b.EdgeColor = 'none';
b.FaceAlpha = .7;

hold on;
bar(0, 'facecolor',[134 156 211]/255, 'edgecolor','none','facealpha', .7);
% bar(0, 'facecolor',[178 198 65]/255, 'edgecolor','none','facealpha', .7);
l = legend('Baseline','Stim','PostStim');

% p.FaceAlpha = 0.5;
xlabel('Recording');
ylabel('Spiking frequency (Hz)');
pbaspect([2,1,1]);
box off;
% p.CData = color_data;
set(gcf, 'units', 'normalized', 'outerposition', [0 0 .7 .7],'color','w');
axis padded
title('210610 Const1 Condition4.5')
exportgraphics(gcf,'/Users/jeremi/Desktop/210514_Const1_Cond5_spcounts.png','resolution',300);