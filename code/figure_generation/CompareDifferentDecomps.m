function CompareDifferentDecomps()

main_decomp = load('./data/processed/main_decomp.mat');

load('./data/processed/alt_decomps.mat','AltEmbedding')

DecompType = {'Diffusion embedding','PCA on connectivity and Phillips genes','PCA on connectivity data only','PCA on gene-expression data only','PCA on first 10 connectvity and gene PCs','Diffusion embedding on mean affinity matrix'};

DecompTypeScatter = {{'','Diffusion embedding: PC1'},{'PCA on connectivity and','Phillips genes: PC1'},{'PCA on connectivity','data only: PC1'},{'PCA on gene-expression','data only: PC1'},{'PCA on first 10 connectvity','and gene PCs: PC1'},{'Diffusion embedding on mean','affinity matrix: PC1'}};

figure('Position',[0 0 1920 1080])

yslice_range = fliplr(218-(87:127))-5;
xslice_range = fliplr(182-(70:90));
plotlabel={'a','b','c','d','e','f','g','h','i','j','k','l'};
 for i = 1:6
     subplot(2,3,i); 
     scatter(main_decomp.pcs_thal(:,1),zscore(AltEmbedding{i}(:,1)),'filled');    
 xlabel('PCA on connectivity and gene data: PC1')
 ylabel(DecompTypeScatter{i})
set(gca,'FontSize',16)
addPlotLabel(plotlabel{i},gca,24,[-.01,0.02])
 end
 
exportgraphics(gcf,'./figure_outputs/S4.png','Resolution',300)
 
 
cmap = turbo(256);

close all

for i = 1:6
 PlotThalGradientSlices(zscore(AltEmbedding{i}(:,1)),main_decomp.used_seed_voxel_coords,turbo(256),[DecompType{i},': PC1 score'],2.1)
  print(['./figure_outputs/Decomp',num2str(i),'.png'],'-dpng','-r300')
end
close all