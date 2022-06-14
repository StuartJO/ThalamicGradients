

MouseGenesPC1 = MouseGeneTable.PC1;

MouseIN = ismember(MouseGenesUpper,GeneIDs_all_valid);

[~,homologue_gene_ind] = ismember(MouseGenesUpper(MouseIN),GeneIDs_all_valid);
MouseGenesUpper(MouseIN);
GeneIDs_all_valid(homologue_gene_ind);

[~,MouseTopGenes] = sort((MouseGenesPC1),'ascend');

idx = (ismember(GeneIDs_all_valid, MouseGenesUpper));

Violin(MouseGenesPC1(MouseIN),1);
hold on
Violin(MouseGenesPC1(~MouseIN),2); 
 
gene_pc1 = coeff{3}(251:end,1);
turbo_cmap = turbo(256);
scatter(MouseGenesPC1(MouseIN),gene_pc1(homologue_gene_ind),60,'filled')
xlabel('Mouse gene PC1 loading')
ylabel({'Human homologue','gene PC1 loading'})
set(gca,'FontSize',24)
%print('MouseHumanPC1Corr','-dpng','-r300')
exportgraphics(gcf,'MouseHumanPC1Corr.png','Resolution',300)

Violin(gene_pc1((ismember(GeneIDs_all_valid, MouseGenesUpper))),1);
hold on
Violin(gene_pc1(~(ismember(GeneIDs_all_valid, MouseGenesUpper))),2); 
xticks(1:2)
xticklabels({'Genes in mouse data','Genes not in mouse data'})
ylabel('PC1 loading')