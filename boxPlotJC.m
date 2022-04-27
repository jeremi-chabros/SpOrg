close all; clc;
plt = 1;

var2plot = gramian_det;
    tl = 'both';
remove_outliers = 0;
caption = "V("+char(949)+")";
% caption = "Modal controllability";
bar_width = 0.2;
jitWidth = bar_width;
sep = 3/2*bar_width;
% gtypes = fieldnames(var2plot);
gtypes = {'WT','HE','KO'};
% gtypes = {'WT','KO'};
pos = 1:size(var2plot.(gtypes{1}));
pos = 3/2*pos;
pos = [pos-sep; pos; pos+sep];
colors = linspecer(3);
colors = [[0, 113, 188]/256; [182, 50, 28]/256; 0.3 0.3 0.3];
ages = {'14','21','28','35'};

mrks = {'o', '^', 's'};
pstats=[];
dzio=[];
dziodzio=[];
sstars=[];

if remove_outliers
    for i = 1:numel(gtypes)
        
        B = var2plot.(gtypes{i});
        X = rmoutliers(B','quartiles');
        X = X';
        var2plot.(gtypes{i}) = X;
        
    end
end

for i = 1:numel(gtypes)
    boxplot(var2plot.(gtypes{i})',...
        'positions', pos(i,:),...
        'widths', bar_width,...
        'colors', colors(i,:),...
        'symbol','');
    set(findobj(gca,'type','line'),'linew',1)
    hold on
    for j = 1:length(pos)
        jitter = linspace(0, bar_width, length(var2plot.(gtypes{i})(j,:)));
        s = scatter(pos(i,j)-jitWidth/2+jitter, (var2plot.(gtypes{i})(j,:)));
        s.Marker = mrks{i};
        s.MarkerFaceColor = colors(i,:);
        s.MarkerEdgeColor = colors(i,:);
        %         s.MarkerEdgeColor = 'k';
%                 s.SizeData = 30;
        s.MarkerFaceAlpha = 0.7;
    end
end

hold on

% Stats testing
for i = 1:size(pos,2)
    
    if length(gtypes)==3
        WT = var2plot.(gtypes{1})(i,:);
        HE = var2plot.(gtypes{2})(i,:);
        KO = var2plot.(gtypes{3})(i,:);
        combs = nchoosek(1:3,2);
        
        
        [pstats(i,2),~,ostats(i,2)] = ranksum(WT, KO, 'tail', tl);
%                 stats_mes = mes(WT',KO','U3');
%                 Fsize(i,1) = stats_mes.U3;
        stats_mes = (ostats(i,2).ranksum)/(length(KO)*length(WT));
        Fsize(i,2) = stats_mes;
        
        [pstats(i,3),~,ostats(i,3)]= ranksum(KO, HE, 'tail', tl);
%                 stats_mes = mes(HE',KO','U3');
%                 Fsize(i,2) =  stats_mes.U3;
        stats_mes = (ostats(i,3).ranksum)/(length(KO)*length(HE));
        Fsize(i,3) = stats_mes;
%         
        [pstats(i,1),~,ostats(i,1)]= ranksum(HE, WT, 'tail', tl);
%                 stats_mes = mes(WT',HE','U3');
%                 Fsize(i,3) = stats_mes.U3;
        stats_mes = (ostats(i,1).ranksum)/(length(HE)*length(WT));
        Fsize(i,1) = stats_mes;


    else
        KO = var2plot.(gtypes{2})(i,:);
        WT = var2plot.(gtypes{1})(i,:);
        combs = nchoosek(1:2,2);
        
        p = ranksum(KO, WT, 'tail', tl);
        pstats(i) = p;
        stats_mes = mes(WT',KO','U3');
        Fsize(i) = stats_mes.U3;
    end
end

for i = 1:4
    for j = 1:length(nchoosek(1:length(gtypes),2))
        if length(gtypes)>2
            sstars{i,j} = [pos(combs(j,1),i),pos(combs(j,2),i)];
        else
            sstars{i} = [pos(1,i), pos(2,i)];
        end
    end
end

G = array2table(Fsize);
G.Properties.VariableNames = {'HE_WT','KO_WT','KO_HE'};
G.Properties.RowNames = {'DIV14','DIV21','DIV28','DIV35'};

S = array2table(pstats);
S.Properties.VariableNames = {'HE_WT','KO_WT','KO_HE'};
S.Properties.RowNames = {'DIV14','DIV21','DIV28','DIV35'};

dzio = {sstars{:}};
dziodzio = pstats(:)';

sigstar(dzio, dziodzio, 0);
hold on
box off;
xticks(pos(2,:))
xticklabels(ages)
xlabel('Age (DIV)')
ylabel(caption)
set(gcf, 'color', 'w');


% Phantom plots for the legend
h1=plot(nan, 'o-', 'color', colors(1,:), 'linewidth', 1,...
    'markerFaceColor', colors(1,:));
hold on
h2=plot(nan, '^-',  'color', colors(2,:), 'linewidth', 1,...
    'markerfacecolor', colors(2,:));
hold on
h3=plot(nan, 's-',  'color', colors(3,:), 'linewidth', 1,...
    'markerfacecolor', colors(3,:));

% Legend
l = legend([h1 h2 h3]);
l.Title.String = 'Genotype:';
l.String = {'WT','HET','KO'};
l.Location = 'northeastoutside';
l.Box = 'off';
ax=gca;
ax.LineWidth = 1;

axis padded

set(gcf,'position',[4 4 10 5]);


if plt
%     exportgraphics(gcf, [caption '.png'], 'resolution', 500)
%     matlab2tikz('erank.tex');
    matlab2tikz('ctrb_ellipsoid.tex','width','5.75in','height','3in');
end
