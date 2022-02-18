for g = 1

PC_gene_coeffs = coeff{3}(251:end,g);

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


plot(EnrichmentLvls,enrichment,'Color',turbo(1))
hold on
scatter(EnrichmentLvls(isSig),enrichment(isSig),30,'filled','MarkerEdgeColor',turbo(1),'MarkerFaceColor',turbo(1))
scatter(EnrichmentLvls(~isSig),enrichment(~isSig),30,'MarkerEdgeColor',turbo(1),'MarkerFaceColor',[1 1 1])

ylabel('Enrichment')
xlabel('Top {\itN} genes')
title(['PC',num2str(g)])
set(gca,'FontSize',24)
ylim([0 10])

end