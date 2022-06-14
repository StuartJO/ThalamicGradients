





var_explained = explained{3};

plot(var_explained,'Color','k','LineWidth',2)
hold on
scatter(1:length(var_explained),var_explained,50,'filled','MarkerEdgeColor','k','MarkerFaceColor','k')
xlim([.5 length(var_explained)])
xlabel('Principal component')
ylabel('% variance explained')
set(gca,'FontSize',24)

spacing = linspace(length(var_explained),6,60);

for i = 1:60
    xlim([.5 spacing(i)+.5])
    print(['./GIF/VAR_EXPLAINED_',num2str(i),'.png'],'-dpng')
end


pc_thal = zscore(score{3}(:,1));
 imagesc(zscore(score{3}))
 
 caxis([min(pc_thal) max(pc_thal)])
 
 colormap(turbo(256))
 
 axis image
 axis off
 
 for i = 1:60
    xlim([.5 spacing(i)+.5])
    %print(['./GIF/SCOREMAT_zoom_',num2str(i),'.png'],'-dpng')
    exportgraphics(gcf,['./GIF/SCOREMAT_zoom_',num2str(i),'.png'])
 end

 
 
 
 
 pc_thal = zscore(coeff{3}(:,1));
 imagesc(zscore(coeff{3}))
 
 caxis([min(pc_thal) max(pc_thal)])
 
 colormap(turbo(256))
 
 axis image
 axis off
  for i = 1:60
    xlim([.5 spacing(i)+.5])
    %print(['./GIF/SCOREMAT_zoom_',num2str(i),'.png'],'-dpng')
    exportgraphics(gcf,['./GIF/COEFFMAT_zoom_',num2str(i),'.png'])
  end
 
   imagesc(zscore(coeff{3}))
   
    caxis([min(pc_thal) max(pc_thal)])
 
 colormap(turbo(256))
 
 exportgraphics(gcf,['./GIF/COEFFMAT_zoom_',num2str(i),'.png'])