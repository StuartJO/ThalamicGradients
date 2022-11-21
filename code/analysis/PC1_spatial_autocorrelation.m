function PC1_spatial_autocorrelation()

load('main_decomp.mat','score','coeff','GeneIDs_valid')
load('TractGeneNorm.mat','GeneData_norm')

pc1_thal = zscore(score(:,1));

surrogates_subsample = csvread('thal_surrogates_NEW.csv')';

Ngenes = size(GeneData_norm,2);

gene_corr_surr_pc1 = zeros(size(surrogates_subsample,2),Ngenes);

for i = 1:Ngenes
gene_corr_surr_pc1(:,i) = corr(GeneData_norm(:,i),surrogates_subsample, ...
            'rows','pairwise','type','pearson');
end

gene_corr_pc1 = corr(pc1_thal,GeneData_norm, ...
            'rows','pairwise','type','pearson');
        
prctile_rank_genes = zeros(Ngenes,1);
        
for i = 1:Ngenes
prctile_rank_genes(i) = mean(gene_corr_surr_pc1(:,i) > gene_corr_pc1(i)); 
end

%significant = prctile_rank_genes < 0.025 | prctile_rank_genes >= 0.975;

p_perm_thal = zeros(Ngenes,1);

for i = 1:Ngenes
    
if gene_corr_pc1(i) > 0
    p_perm_thal(i) = mean(gene_corr_surr_pc1(:,i) > gene_corr_pc1(i));
else
    p_perm_thal(i) = mean(gene_corr_surr_pc1(:,i) < gene_corr_pc1(i));
end

end

significant = p_perm_thal < .05;

PC_gene_coeffs = coeff(251:end,1);

PC_gene_coeffs_significant = PC_gene_coeffs(significant);
GeneIDs_significant = GeneIDs_valid(significant);

[PC_genes_positive,PC_gene_coeffs_sorted_ind] = sort(PC_gene_coeffs_significant,'descend');

MostPositive = GeneIDs_significant(PC_gene_coeffs_sorted_ind(1:100));
[PC_genes_negative,PC_gene_coeffs_sorted_negative_ind] = sort(PC_gene_coeffs_significant,'ascend');
MostNegative = GeneIDs_significant(PC_gene_coeffs_sorted_negative_ind(1:100));

writecell([MostPositive num2cell(PC_genes_positive(1:100))],'./data/processed/HumanMostPositiveSpinTested.csv')
writecell([MostNegative num2cell(PC_genes_negative(1:100))],'./data/processed/HumanMostNegativeSpinTested.csv')

writecell([GeneIDs_valid num2cell(PC_gene_coeffs) num2cell(p_perm_thal)],'./data/processed/AllHumanSpinTested.csv')