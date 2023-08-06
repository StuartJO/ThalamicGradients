

load('main_decomp.mat','GeneNames_human')

for i = 1:3
surrogates = csvread(['thal_surrogates_PC',num2str(i),'.csv']);
AutoCorrReslt = PC_spatial_autocorrelation(score(:,i),coeff(251:end,i),surrogates');

writecell([AutoCorrReslt.MostPositive num2cell(AutoCorrReslt.PC_genes_positive(1:100))],['./data/processed/HumanMostPositiveSpinTested_PC',num2str(i),'.csv'])
writecell([AutoCorrReslt.MostNegative num2cell(AutoCorrReslt.PC_genes_negative(1:100))],['./data/processed/HumanMostNegativeSpinTested_PC',num2str(i),'.csv'])
writecell([GeneNames_human num2cell(coeff(251:end,i)) num2cell(AutoCorrReslt.p_perm)],['./data/processed/AllHumanSpinTested_PC',num2str(i),'.csv'])

end

