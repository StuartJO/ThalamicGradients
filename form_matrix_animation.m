data = data_type{1};
d = squareform(pdist(data));
seed_linkage = linkage(data,'ward');
seed_order = optimalleaforder(seed_linkage,d);
tract_data = data_type{1};
d = squareform(pdist(tract_data'));
tract_linkage = linkage(tract_data','ward');
tract_order = optimalleaforder(tract_linkage,d);
gene_data = data_type{2};
d = squareform(pdist(gene_data'));
gene_linkage = linkage(gene_data','ward');
gene_order = optimalleaforder(gene_linkage,d);



imagesc(tract_data(seed_order,tract_order))
xticks([])
yticks([])
colormap(parula(1000))

order_tract_data = tract_data(seed_order,tract_order);
order_tract_data_orig=order_tract_data;


order_tract_data(2:end,:) = NaN;
tract_plot = pcolor(order_tract_data);
tract_plot.EdgeAlpha=0;
set(gca,'YDir','reverse')
xlim([1 10])
ylim([-3 6])

xlinspace = [10 10 linspace(10,251,118)];

ylinspace = [xlinspace(1:30)-4 linspace(xlinspace(31)-4,921,90)];

ylinspace2 = linspace(-3,1,120);

rowspace = [2:30 round(linspace(31,921,91))];
% yaxisscale = (1:90).^10;
% ylinspace = rescale(yaxisscale./max(yaxisscale),6,921);
axis square
axis off
for i = 1:120
    tract_plot.CData(1:rowspace(i),:) = order_tract_data_orig(1:rowspace(i),:);
    xlim([1 xlinspace(i)])
    ylim([ylinspace2(i) ylinspace(i)])
    pause(.1)
end






 pc_thal = zscore(score{3}(:,1));
 
 imagesc(zscore(score{3}))
 
 caxis([min(pc_thal) max(pc_thal)])
 
 colormap(turbo(256))
 
  imagesc(sort(pc_thal))
 colormap(turbo(256))
 
 pc_cort_gene = zscore(coeff{3}(:,1));
 
 
    imagesc(zscore(coeff{3}))
 colormap(turbo(256))
  caxis([min(pc_cort_gene) max(pc_cort_gene)])
  
  sum(ThalSeedAvg(:,1:250)>0,2);
  
  sum(ThalSeedAvg(:,1:250),2);
  
 pc_thal = zscore(score{3}(:,1));  
imagesc(score{3}(1:20,1:5))
colormap(turbo(256))


pc_cort_loadings = zscore(coeff{3}(1:250,1));

imagesc(pc_cort_loadings(1:250))
colormap(turbo(256))
axis off
exportgraphics(gcf,'HIGHRES_PC1_cort_loadings.png','Resolution',300)

pc_gene_loadings = zscore(coeff{3}(251:end,1));

imagesc(pc_gene_loadings)
colormap(turbo(256))
axis off
exportgraphics(gcf,'HIGHRES_PC1_gene_loadings.png','Resolution',300)


pc_thal = zscore(score{3}(:,1));  
zscore_scores = zscore(score{3});
imagesc(zscore_scores(randi(921,1,20),1:20))
colormap(turbo(256))
caxis([min(pc_thal) max(pc_thal)])
axis square
yticks(1:20)
xticks(1:20)
set(gca,'XAxisLocation', 'top','TickLength',[0 0])
for i = 1:20
   x_ticklabels{i} = ['PC',num2str(i)]; 
   y_ticklabels{i} = ['Seed ',num2str(i)]; 
end
yticklabels(y_ticklabels)
xticklabels(x_ticklabels)
xtickangle(45)
set(gca,'FontSize',20)
print('./ScoreMatZoom.png','-dpng','-r300')



 pc_cort_gene = zscore(coeff{3});
 
pc_thal = zscore(score{3}(:,1));  
zscore_scores = zscore(score{3});
%imagesc(zscore_scores(randi(921,1,20),1:20))
imagesc(pc_cort_gene(1:300,1:100))
colormap(turbo(256))
caxis([min(pc_cort_gene(:,1)) max(pc_cort_gene(:,1))])
axis image
yticks(1:300)
xticks(1:100)
set(gca,'XAxisLocation', 'top','TickLength',[0 0])
for i = 1:300
  x_ticklabels{i} = ['PC',num2str(i)]; 
   if i <=250
   y_ticklabels{i} = ['ROI ',num2str(i)]; 
   else
     y_ticklabels{i} = ['Gene ',num2str(i-250)];   
   end
end
yticklabels(y_ticklabels)
xticklabels(x_ticklabels(1:100))
xtickangle(45)
set(gca,'FontSize',2)
print('./CoeffMatZoom.png','-dpng','-r1200')


