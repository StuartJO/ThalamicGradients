load('./data/processed/decomp_rand500.mat')

data = pca_input;
d = squareform(pdist(data));
seed_linkage = linkage(data,'ward');
seed_order = optimalleaforder(seed_linkage,d);
tract_data = pca_input(:,1:250);
d = squareform(pdist(tract_data'));
tract_linkage = linkage(tract_data','ward');
tract_order = optimalleaforder(tract_linkage,d);
gene_data = pca_input(:,251:end);
d = squareform(pdist(gene_data'));
gene_linkage = linkage(gene_data','ward');
gene_order = optimalleaforder(gene_linkage,d);
imagesc(tract_data(seed_order,tract_order))
xticks([])
yticks([])
colormap(parula(1000))
colorbar('Location','SouthOutside')
set(gca,'FontSize',24)
print('./figure_outputs/ConnMatrix.svg','-dsvg')

figure

imagesc(gene_data(seed_order,gene_order))
colormap(viridis)
xticks([])
yticks([])
colorbar('Location','SouthOutside')
set(gca,'FontSize',24)

print('./figure_outputs/GeneMatrix.svg','-dsvg')

figure

bar_cmap = brewermap(512,'RdBu');
pos_bar_cmap = flipud(bar_cmap(1:256,:));

joint_data = [tract_data(seed_order,tract_order) gene_data(seed_order,gene_order)];
imagesc(joint_data)
colormap(pos_bar_cmap)
xticks([])
yticks([])
colorbar('Location','SouthOutside')
set(gca,'FontSize',24)
print('./figure_outputs/JointMatrix.svg','-dsvg')

figure('Position',[349 98 180 1068])
imagesc(zscore(score(seed_order,1)))
colormap(turbo(256))
xticks([])
yticks([])

print('./figure_outputs/PC1ScoreMatrix.svg','-dsvg')

figure('Position',[349 98 180 1068])
imagesc(zscore(coeff([tract_order gene_order+length(tract_order)],1)))
colormap(turbo(256))
xticks([])
yticks([])
print('./figure_outputs/PC1CoeffMatrix.svg','-dsvg')