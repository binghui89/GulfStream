load box_plot;
share = [share_matrix_single; share_matrix_pa];
xticks_matrix_single = repmat(11: 19, 19188, 1);
xticks_matrix_pa = repmat(21: 29, 20, 1);
xticks = [xticks_matrix_single; xticks_matrix_pa];

positions = 0.3333: 8.3333;
positions = [positions; positions + 0.3333];

bh = boxplot(share(:), int2str(xticks(:)), 'positions', positions(:),...
    'width', 0.3);
for i=1:size(bh,2) % <- # graphics handles/x
     set(bh(:,i),'linewidth',1);
%      disp(sprintf('working on component: %3d = %s',i,get(bh(i,1),'tag')));
%      pause(.5);
end


hold on;
for x = 0: 9
   line([x, x], [-1, 1], 'color', [.5, .5, .5], 'LineStyle', '--'); 
end
xlim([0, 9]);
xlabel_positions = mean(positions, 1);
set(gca,'xtick',xlabel_positions);
set(gca,'xticklabel',...
    {'GEN', 'TL', 'MOOR', 'DPLY', 'DEV', 'SUB', 'CON', 'FOM', 'INS'});
color = repmat(['k', 'w'], 1, 9);
h = findobj(gca,'Tag','Box');
for j=1:length(h)
   patch(get(h(j),'XData'),get(h(j),'YData'),color(j),'FaceAlpha',.5);
end

c = get(gca, 'Children');
hleg1 = legend(c(1:2), 'Single grid', 'Portfolio');

set(gca, 'FontSize', 14);
hold off;