function AutoCorrReslt = PC_spatial_autocorrelation(pc_score,pc_gene_coeffs,TractStruct,surrogates)

%load('main_decomp.mat','GeneNames_human')
%load('TractGeneNorm.mat','GeneData_norm','GeneNames_human')

%pc = zscore(score(:,1));

%surrogates = csvread('thal_surrogates_NEW.csv')';

Ngenes = size(TractStruct.GeneData_norm,2);

gene_corr_surr_pc = zeros(size(surrogates,2),Ngenes);

for i = 1:Ngenes
gene_corr_surr_pc(:,i) = corr(TractStruct.GeneData_norm(:,i),surrogates, ...
            'rows','pairwise','type','pearson');
end

gene_corr_pc = corr(pc_score,TractStruct.GeneData_norm, ...
            'rows','pairwise','type','pearson');
        
prctile_rank_genes = zeros(Ngenes,1);
        
for i = 1:Ngenes
prctile_rank_genes(i) = mean(gene_corr_surr_pc(:,i) > gene_corr_pc(i)); 
end

%significant = prctile_rank_genes < 0.025 | prctile_rank_genes >= 0.975;

p_perm = zeros(Ngenes,1);

for i = 1:Ngenes
    
if gene_corr_pc(i) > 0
    p_perm(i) = mean(gene_corr_surr_pc(:,i) > gene_corr_pc(i));
else
    p_perm(i) = mean(gene_corr_surr_pc(:,i) < gene_corr_pc(i));
end

end

significant = p_perm < .05;

AutoCorrReslt.p_perm = p_perm;

PC_gene_coeffs_significant = pc_gene_coeffs(significant);
GeneIDs_significant = TractStruct.GeneNames_human(significant);

[PC_genes_positive,PC_gene_coeffs_sorted_ind] = sort(PC_gene_coeffs_significant,'descend');

MostPositive = GeneIDs_significant(PC_gene_coeffs_sorted_ind(1:100));
[PC_genes_negative,PC_gene_coeffs_sorted_negative_ind] = sort(PC_gene_coeffs_significant,'ascend');
MostNegative = GeneIDs_significant(PC_gene_coeffs_sorted_negative_ind(1:100));

AutoCorrReslt.PC_genes_positive = PC_genes_positive;
AutoCorrReslt.MostPositive = MostPositive;

AutoCorrReslt.PC_genes_negative = PC_genes_negative;
AutoCorrReslt.MostNegative = MostNegative;

%writecell([MostPositive num2cell(PC_genes_positive(1:100))],'./data/processed/HumanMostPositiveSpinTested.csv')
%writecell([MostNegative num2cell(PC_genes_negative(1:100))],'./data/processed/HumanMostNegativeSpinTested.csv')
%writecell([GeneNames_human num2cell(pc_gene_coeffs) num2cell(p_perm_thal)],'./data/processed/AllHumanSpinTested.csv')