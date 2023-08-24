%% Get significantly correlated genes

addpath(genpath('./'))

main_decomp = load('main_decomp.mat');
main_genedata = load('TractGeneNorm_rand500.mat');
for i = 1:3
surrogates = csvread(['thal_surrogates_PC',num2str(i),'.csv']);  
pc_gene_coeffs = main_decomp.coeff(251:end,i);
AutoCorrReslt = PC_spatial_autocorrelation(main_decomp.score(:,i),pc_gene_coeffs,main_genedata,surrogates');

writecell([AutoCorrReslt.MostPositive num2cell(AutoCorrReslt.PC_genes_positive(1:100))],['./data/processed/PC',num2str(i),'_HumanMostPositiveSpinTested.csv'])
writecell([AutoCorrReslt.MostNegative num2cell(AutoCorrReslt.PC_genes_negative(1:100))],['./data/processed/PC',num2str(i),'_HumanMostNegativeSpinTested.csv'])
writecell([main_decomp.GeneNames_human num2cell(pc_gene_coeffs) num2cell(AutoCorrReslt.p_perm)],['./data/processed/PC',num2str(i),'_AllHumanSpinTested.csv'])

writecell(AutoCorrReslt.MostPositive,['./data/processed/PC',num2str(i),'_HumanMostPositiveSpinTested.txt'])
writecell(AutoCorrReslt.MostNegative,['./data/processed/PC',num2str(i),'_HumanMostNegativeSpinTested.txt'])

end