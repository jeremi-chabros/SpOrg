var2plot = modal_ctrb;
plt = 1;
g = 'KO';
clr = colors(3,:);
cap = 'Modal controllability';
for i = 15:-1:1
    lh = plot(var2plot.(g)(:,i), 's-', 'color', [.5 .5 .5]);
    lh.Color = [lh.Color .5];
    lh.MarkerFaceColor = clr;
    hold on;
end

hold on
ebar = errorbar([1 2 3 4], median(var2plot.(g)'), iqr(var2plot.(g),2));
ebar.Color = clr;
ebar.LineWidth = 1.5;
box off;
xticks([1 2 3 4])
xticklabels({'14' '21' '28' '35'});
axis padded
set(gcf,'position',[4 4 10 5]);
ylabel(cap)
xlabel('Age (DIV)')
ax=gca;
ax.LineWidth = 1;

if plt
   
    exportgraphics(gcf, [cap '_' g '.png'], 'resolution', 500)
    
end