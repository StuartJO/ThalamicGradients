function PlotMousePC1vsXcoord()

load('MouseThalROICoord.mat','MouseThalROICoord')

figure('Position',[168   187   823   684])
scatter(MouseThalROICoord,thal_pc,80,mouse_region_PC1_color,'filled')
hold on
for i = 1:length(mouse_medpos)
    text(MouseThalROICoord((i)),thal_pc((i))+.1,structInfo.acronym{(i)+87},'FontSize',12,'HorizontalAlignment','center');
end
set(gca,'FontSize',20)
xlabel('CCFv3 {\itx}-coordinate')
ylabel('mPC1 score')
%addPlotLabel('a')
print(['./figure_outputs/ThalNucleiNamePos.svg'],'-dsvg')
print(['./figure_outputs/ThalNucleiNamePos.png'],'-dpng')