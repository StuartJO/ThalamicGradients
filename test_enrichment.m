MouseGenesPC1 = MouseGeneTable.PC1;

idx = find(ismember(GeneIDs_all_valid, MouseGenesUpper));

[~,MouseTopGenes] = sort((MouseGenesPC1),'descend');

idx = find(ismember(GeneIDs_all_valid, MouseGenesUpper(MouseTopGenes(1:100))));

turbo_cmap = turbo(2);

for grad = 1

PC_gene_coeffs = coeff{3}(251:end,grad);

[~,PC_gene_coeffs_sorted] = sort((PC_gene_coeffs),'descend');

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

%figure('Position',[238 304 758 413])
plot(EnrichmentLvls,enrichment,'Color',turbo_cmap(2,:))
hold on
scatter(EnrichmentLvls(isSig),enrichment(isSig),30,'filled','MarkerEdgeColor',turbo_cmap(2,:),'MarkerFaceColor',turbo_cmap(2,:))
scatter(EnrichmentLvls(~isSig),enrichment(~isSig),30,'MarkerEdgeColor',turbo_cmap(2,:),'MarkerFaceColor',[1 1 1])

ylabel('Enrichment')
xlabel({'Top {\itN} genes with ','largest positive PC1 loadings'})
%title(['PC',num2str(g)])
set(gca,'FontSize',24)
%ylim([0 10])
end
exportgraphics(gcf,['./enrichment_positive.png'],'Resolution',300)

clf

[~,MouseTopGenes] = sort((MouseGenesPC1),'ascend');

idx = find(ismember(GeneIDs_all_valid, MouseGenesUpper(MouseTopGenes(1:100))));

turbo_cmap = turbo(2);

for grad = 1

PC_gene_coeffs = coeff{3}(251:end,grad);

[~,PC_gene_coeffs_sorted] = sort((PC_gene_coeffs),'ascend');

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

%figure('Position',[238 304 758 413])
plot(EnrichmentLvls,enrichment,'Color',turbo_cmap(1,:))
hold on
scatter(EnrichmentLvls(isSig),enrichment(isSig),30,'filled','MarkerEdgeColor',turbo_cmap(1,:),'MarkerFaceColor',turbo_cmap(1,:))
scatter(EnrichmentLvls(~isSig),enrichment(~isSig),30,'MarkerEdgeColor',turbo_cmap(1,:),'MarkerFaceColor',[1 1 1])

ylabel('Enrichment')
xlabel({'Top {\itN} genes with ','largest negative PC1 loadings'})
%title(['PC',num2str(g)])
set(gca,'FontSize',24)
%ylim([0 10])
end
exportgraphics(gcf,['./enrichment_negative.png'],'Resolution',300)


gene_loadings = coeff{3}(251:end,grad);

[gene_loadings_order,order_max] = sort(gene_loadings,'descend');

gene_loadings_nan = [gene_loadings gene_loadings];

gene_loadings_nan(order_max(101:end),:)=NaN;

h = imagesc2(gene_loadings_nan);

colormap(turbo(256))

caxis([min(gene_loadings) max(gene_loadings)])
axis off 
axis normal
print(['./GENES_top.png'],'-dpng')

[~,order_min] = sort(gene_loadings,'ascend');

gene_loadings_nan = [gene_loadings gene_loadings];

gene_loadings_nan(order_min(101:end),:)=NaN;
h = imagesc2(gene_loadings_nan);
colormap(turbo(256))

caxis([min(gene_loadings) max(gene_loadings)])
axis off 
axis normal
print(['./GENES_bottom.png'],'-dpng')

gene_loadings_nan = [gene_loadings gene_loadings];
h = imagesc2(gene_loadings_nan);
colormap(turbo(256))

caxis([min(gene_loadings) max(gene_loadings)])
axis off 
axis normal
print(['./GENES_all.png'],'-dpng')