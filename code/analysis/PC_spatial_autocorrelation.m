function AutoCorrReslt = PC_spatial_autocorrelation(pc_score, pc_gene_coeffs, TractStruct, surrogates)

% This function calculates spatial autocorrelation of principal components (PCs)
% with respect to gene expression data and identifies significant correlations.

% Inputs:
%   pc_score: Principal component scores matrix. Each row corresponds to a sample, 
%             and each column corresponds to a principal component score.
%   pc_gene_coeffs: Coefficients matrix representing the correlation between
%                   principal components and gene expression data. Each row
%                   corresponds to a gene, and each column corresponds to a principal component.
%   TractStruct: A structure containing gene expression data.
%                TractStruct.GeneData_norm: Normalized gene expression data matrix.
%                TractStruct.GeneNames_human: Cell array of gene names.
%   surrogates: Matrix of surrogate gene expression data for statistical comparison.

% Outputs:
%   AutoCorrReslt: A structure containing the following fields:
%       - p_perm: Array of p-values obtained from permutation tests for each gene.
%       - PC_genes_positive: Cell array of gene names with the highest positive correlation
%                           coefficients with the principal components.
%       - MostPositive: List of gene names with the highest positive correlation coefficients.
%       - PC_genes_negative: Cell array of gene names with the highest negative correlation
%                           coefficients with the principal components.
%       - MostNegative: List of gene names with the highest negative correlation coefficients.

Ngenes = size(TractStruct.GeneData_norm, 2);

gene_corr_surr_pc = zeros(size(surrogates, 2), Ngenes);

% Calculate gene correlation with surrogate gene expression data
for i = 1:Ngenes
    gene_corr_surr_pc(:, i) = corr(TractStruct.GeneData_norm(:, i), surrogates, ...
        'rows', 'pairwise', 'type', 'pearson');
end

% Calculate gene correlation with actual principal component scores
gene_corr_pc = corr(pc_score, TractStruct.GeneData_norm, ...
    'rows', 'pairwise', 'type', 'pearson');

prctile_rank_genes = zeros(Ngenes, 1);

% Calculate the percentile rank of gene correlations using surrogate data
for i = 1:Ngenes
    prctile_rank_genes(i) = mean(gene_corr_surr_pc(:, i) > gene_corr_pc(i)); 
end

p_perm = zeros(Ngenes, 1);

% Calculate permutation test p-values
for i = 1:Ngenes
    if gene_corr_pc(i) > 0
        p_perm(i) = mean(gene_corr_surr_pc(:, i) > gene_corr_pc(i));
    else
        p_perm(i) = mean(gene_corr_surr_pc(:, i) < gene_corr_pc(i));
    end
end

% Identify significant correlations based on a threshold of 0.05
significant = p_perm < 0.05;

% Store results in the output structure
AutoCorrReslt.p_perm = p_perm;

% Select genes with significant positive and negative correlations
PC_gene_coeffs_significant = pc_gene_coeffs(significant);
GeneIDs_significant = TractStruct.GeneNames_human(significant);

[PC_genes_positive, PC_gene_coeffs_sorted_ind] = sort(PC_gene_coeffs_significant, 'descend');
MostPositive = GeneIDs_significant(PC_gene_coeffs_sorted_ind(1:100));

[PC_genes_negative, PC_gene_coeffs_sorted_negative_ind] = sort(PC_gene_coeffs_significant, 'ascend');
MostNegative = GeneIDs_significant(PC_gene_coeffs_sorted_negative_ind(1:100));

AutoCorrReslt.PC_genes_positive = PC_genes_positive;
AutoCorrReslt.MostPositive = MostPositive;

AutoCorrReslt.PC_genes_negative = PC_genes_negative;
AutoCorrReslt.MostNegative = MostNegative;
