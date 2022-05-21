load('Sub76_ThalData.mat')

EntrezIDs_all = dlmread('AHBAEntrez.txt');

GeneIDs_all_table = readtable('AHBAgeneSymbol.txt','ReadVariableNames',false);

GeneIDs_all = GeneIDs_all_table.Var1;

[~,GeneIDs_all_Gryglewski_ind] = ismember(GenesKept,EntrezIDs_all);

GeneIDs_all_valid = GeneIDs_all(GeneIDs_all_Gryglewski_ind);

MouseGeneTable = readtable('MouseThalGenes.xlsx');

MouseGenes = MouseGeneTable.GeneSymbol;

MouseGenesUpper = upper(MouseGenes);

idx = find(ismember(GeneIDs_all_valid, MouseGenesUpper));

MouseGenesValid = GeneIDs_all_valid(idx);

% data_type{1} = zscore(ThalSeedAvg(:,1:250));
% data_type{2} = zscore(SeedGene_kept);

data_type{1} = BF_NormalizeMatrix(ThalSeedAvg(:,1:250),'scaledSigmoid');
data_type{2} = BF_NormalizeMatrix(SeedGene_kept,'scaledSigmoid');

data_type{3} = [data_type{1} data_type{2}];
data_type{4} = data_type{2}(:,idx);
data_type{5} = [data_type{1} data_type{2}(:,idx)];




for PCtype = 1:5
[coeff{PCtype},score{PCtype},~,~,explained{PCtype}] = pca(data_type{PCtype});
end

figure

grad = 1;
PCtype = 3;

for grad = [1 2]
    for PCtype = [1 3]

PCtypenames = {'tract','gene','joint'};
XYZ = 1;

load('MNI_Seed_voxelData.mat','seeds_vox','seed_coords')

seed_voxel_coords = seeds_vox(logical(seed_ind),:);

pc_thal = zscore(score{PCtype}(:,grad));
    
%PlotThalGradient2(pc_thal,seed_voxel_coords,turbo(256),[data_type_name{j},': PCA Gradient ',num2str(g)],2.1)
cmap = turbo(256);
%print(['./NewFigures/PCA_',data_type_savename{j},'_THAL_PC',num2str(g),'.png'],'-dpng','-r300')
PlotThalGradient3(pc_thal,seed_voxel_coords,turbo(256),['Thalamic seed PC',num2str(grad),' score'],2.1)
print(['./NewFigures/',PCtypenames{PCtype},'_Thalamus_PC',num2str(grad),'.png'],'-dpng','-r300')
figure('Position',[162   233   713   592])
seed_mni_coords = seed_coords(logical(seed_ind))';
SPATIAL_DIR = {'Medial-Lateral','Anterior-posterior','Dorsal-ventral'};
      pc_thal = zscore(score{PCtype}(:,grad));
    s = scatter(seed_mni_coords(:,XYZ),pc_thal,'filled','MarkerFaceAlpha',.25,'MarkerFaceColor',cmap(50,:));
    [RHO,pval] = corr(seed_mni_coords(:,XYZ),pc_thal,'Type','Spearman');
    %title([SPATIAL_DIR{i},' ',data_type_name{j},' \rho = ',num2str(RHO)])
    xlabel({'Medial-lateral position','(MNI x-axis coordinate)'})
    ylabel(['Thalamic seed PC',num2str(grad),' score'])
    set(gca,'FontSize',24,'XDir','reverse')
    text(-0.5,1.5,['{\itr_{s}} = ',num2str(round(RHO,3))],'FontSize',24)
    pval_string = num2str(pval,3);
    
    pval_m = pval_string(1:find(pval_string=='e')-1);
    if length(pval_m) < 4
        pval_m = [pval_m,'0'];
    end
    pval_n = pval_string(find(pval_string=='e')+1:end);
    text(-0.5,1,['{\itp} = ',pval_m,'\times10^{',pval_n,'}'],'FontSize',24)
    print(['./NewFigures/',PCtypenames{PCtype},'_','PC',num2str(grad),'_medlat_corr.png'],'-dpng','-r300')
    
    
load('fsaverage_surface_data.mat')
surface.vertices = lh_inflated_verts;
surface.faces = lh_faces;
cmap = turbo(256);

pc = zscore(coeff{PCtype}(1:250,grad));
    
figure('Position',[461   462   560   325])
ax_sub1 = axes('Position',[0.005 .33 .49 .66]);
p = plotSurfaceROIBoundary(surface,lh_rand500,pc,'midpoint',cmap,1,2);
camlight(80,-10);
camlight(-80,-10);
view([-90 0])

axis off
axis image

ax_sub2 = axes('Position',[.505 .33 .489 .66]);
[~,~,~,~,~,orig_data_climits] = plotSurfaceROIBoundary(surface,lh_rand500,pc,'midpoint',cmap,1,2);
camlight(80,-10);
camlight(-80,-10);
view([90 0])
axis off
axis image

c = colorbar('Location','southoutside');
set(c, 'xlim', orig_data_climits,'Position',[.1 .23 .8 .05],'FontSize',20);
c.Label.String = ['Cortical region PC',num2str(grad),' loading'];

print(['./NewFigures/Cortical_',PCtypenames{PCtype},'_','PC',num2str(grad),'.png'],'-dpng','-r300')

end
end

for grad = 1

PC_gene_coeffs = coeff{3}(251:end,grad);

[~,PC_gene_coeffs_sorted] = sort(PC_gene_coeffs,'descend');

hit_list = GeneIDs_all_valid(idx);

full_gene_list = GeneIDs_all_valid;

EnrichmentLvls = 1:100;

enrichment = zeros(1,length(EnrichmentLvls));
p = enrichment;

for i = 1:length(EnrichmentLvls)

top_genes = GeneIDs_all_valid(PC_gene_coeffs_sorted(1:EnrichmentLvls(i)));
%top_genes = MouseGenesValid(PC_gene_coeffs_sorted(1:i));

[enrichment(i),p(i)] = calculate_enrichment(hit_list, top_genes, full_gene_list);

end

[isSig,sigLevel,pvals_corr] = BF_FDR(p,.05,1);
% subplot(1,3,g)
% 
% yyaxis left
% plot(enrichment)
% ylabel('Enrichment')
% 
% yyaxis right 
% plot(p)
% %set(gca, 'YScale', 'log')
% ylabel('p-value')

figure('Position',[238 304 758 413])
plot(EnrichmentLvls,enrichment,'Color',turbo(1))
hold on
scatter(EnrichmentLvls(isSig),enrichment(isSig),30,'filled','MarkerEdgeColor',turbo(1),'MarkerFaceColor',turbo(1))
scatter(EnrichmentLvls(~isSig),enrichment(~isSig),30,'MarkerEdgeColor',turbo(1),'MarkerFaceColor',[1 1 1])

ylabel('Enrichment')
xlabel('Top {\itN} genes')
%title(['PC',num2str(g)])
set(gca,'FontSize',24)
%ylim([0 10])
print(['./NewFigures/Enrichment_PC1.png'],'-dpng','-r300')
end




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
exportgraphics(gcf,'ConnMatrix.png','Resolution',300)
colorbar('Location','SouthOutside')
set(gca,'FontSize',24)
exportgraphics(gcf,'ConnMatrixColorbar.png','Resolution',300)

imagesc(gene_data(seed_order,gene_order))
colormap(viridis)
xticks([])
yticks([])
exportgraphics(gcf,'GeneMatrix.png','Resolution',300)
colorbar('Location','SouthOutside')
set(gca,'FontSize',24)
exportgraphics(gcf,'GeneMatrixColorbar.png','Resolution',300)

joint_data = [tract_data(seed_order,tract_order) gene_data(seed_order,gene_order)+2];
imagesc(joint_data)
parula_cmap = parula(1000);
viridis_cmap = viridis(1000);
colormap([parula_cmap;repmat(parula_cmap(1000,:),500,1); repmat(viridis_cmap(1,:),499,1); viridis_cmap])
xticks([])
yticks([])
exportgraphics(gcf,'JointMatrix.png','Resolution',300)



imagesc(zscore(score{3}(seed_order,1)))
colormap(turbo(256))
xticks([])
yticks([])
exportgraphics(gcf,'PC1ScoreMatrix.png','Resolution',300)

imagesc(zscore(coeff{3}([tract_order gene_order+length(tract_order)],1)))
exportgraphics(gcf,'PC1CoeffsMatrix.png','Resolution',300)

figure
addpath 'F:\Documents\NetGenExample'
CMAP = cmocean('Balance');

d = squareform(pdist(C));
c_linkage = linkage(C,'ward');
c_order = optimalleaforder(c_linkage,d);

imagesc(C(c_order,c_order));
caxis([0 1])
colormap(CMAP)
axis square
xticks([])
yticks([])
exportgraphics(gcf,'CorrMat.png','Resolution',600)
colorbar('Location','SouthOutside')
set(gca,'FontSize',24)
exportgraphics(gcf,'CorrMat_cbar.png','Resolution',300)



PlotThalGradient3(zscore(gene_data(:,2)),seed_voxel_coords,viridis(256),['Gene'],2.1)
print(['./NewFigures/Gene2.png'],'-dpng','-r300')

PlotThalGradient3(zscore(gene_data(:,1)),seed_voxel_coords,viridis(256),['Gene'],2.1)
print(['./NewFigures/Gene1.png'],'-dpng','-r300')

PlotThalGradient3(zscore(gene_data(:,3)),seed_voxel_coords,viridis(256),['Gene'],2.1)
print(['./NewFigures/Gene3.png'],'-dpng','-r300')

PlotThalGradient3(zscore(gene_data(:,4)),seed_voxel_coords,viridis(256),['Gene'],2.1)
print(['./NewFigures/Gene4.png'],'-dpng','-r300')