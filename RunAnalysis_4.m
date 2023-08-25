%% Get significantly correlated genes

addpath(genpath('./'))

decomp_rand500 = load('./data/processed/decomp_rand500.mat');
main_genedata = load('./data/processed/TractGeneNorm_rand500.mat');
for i = 1:3
surrogates = csvread(['./data/processed/thal_surrogates_PC',num2str(i),'.csv']);  
pc_gene_coeffs = decomp_rand500.coeff(251:end,i);
AutoCorrReslt = PC_spatial_autocorrelation(decomp_rand500.score(:,i),pc_gene_coeffs,main_genedata,surrogates');

writecell([AutoCorrReslt.MostPositive num2cell(AutoCorrReslt.PC_genes_positive(1:100))],['./data/processed/PC',num2str(i),'_HumanMostPositiveSpinTested.csv'])
writecell([AutoCorrReslt.MostNegative num2cell(AutoCorrReslt.PC_genes_negative(1:100))],['./data/processed/PC',num2str(i),'_HumanMostNegativeSpinTested.csv'])
writecell([decomp_rand500.GeneNames_human num2cell(pc_gene_coeffs) num2cell(AutoCorrReslt.p_perm)],['./data/processed/PC',num2str(i),'_AllHumanSpinTested.csv'])

if i == 1
    % Copy a version to ./data/gene_data/gene_lists
   writecell([decomp_rand500.GeneNames_human num2cell(pc_gene_coeffs) num2cell(AutoCorrReslt.p_perm)],'./data/gene_data/gene_lists/AllHumanSpinTested.csv') 
end

writecell(AutoCorrReslt.MostPositive,['./data/processed/PC',num2str(i),'_HumanMostPositiveSpinTested.txt'])
writecell(AutoCorrReslt.MostNegative,['./data/processed/PC',num2str(i),'_HumanMostNegativeSpinTested.txt'])

end