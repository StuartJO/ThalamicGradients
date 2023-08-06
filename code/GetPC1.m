% CompiledTractData contains the following variables:
% GenesKept = Entrez IDs of the genes obtained from Gryglewski
% SUB = IDs of HCP subjects whose tractography data was used
% SeedGene_kept = matrix of seeds-by-gene-expression
% ThalSeed = matrix of seeds-by-ROI for each SUB. ROIs=500 because left and
% right hemispheres are included
% ThalSeedAvg = ThalSeed averaged (mean) across subjects
% seed_ind = index of which of the 1811 seeds were kept
load('CompiledTractGeneData.mat')

% Load in gene names and IDs for the AHBA data. These lists were obtained
% from Burt et al., 2020 Nature Neuroscience
EntrezIDs_all = dlmread('AHBAEntrez.txt');
GeneIDs_burt_table = readtable('AHBAgeneSymbol.txt','ReadVariableNames',false);

GeneIDs_burt = GeneIDs_burt_table.Var1;

% Find the names of the genes that existed in Gryglewski et al and Burt et
% al
[~,GeneIDs_burt_Gryglewski_ind] = ismember(GenesKept,EntrezIDs_all);

GeneNames = GeneIDs_burt(GeneIDs_burt_Gryglewski_ind);

% Load in the Phillips et al genes

PhillipsGeneTable = readtable('PhillipsMouseThalGenes.xlsx');

GeneIDs_phillips = PhillipsGeneTable.GeneSymbol;

GeneIDs_phillips_upper = upper(GeneIDs_phillips);

PhillipsMouse_HumanGenesOverlap_idx = find(ismember(GeneNames, GeneIDs_phillips_upper));

TractData_norm = BF_NormalizeMatrix(ThalSeedAvg(:,1:250),'scaledSigmoid');
GeneData_norm = BF_NormalizeMatrix(SeedGene_kept,'scaledSigmoid');

TractData_GeneData_norm = [TractData_norm GeneData_norm];

[coeff,score,~,~,explained] = pca(TractData_GeneData_norm);

load('./data/ancillary/MNI_Seed_voxelData.mat','seed_vox_coords','seed_mni_coords')

used_seed_voxel_coords = seed_vox_coords(logical(seed_ind),:);
used_seed_mni_coords = seed_mni_coords(logical(seed_ind),:);

save('./data/processed/TractGeneNorm.mat','TractData_norm','GeneData_norm')

pc1_thal = zscore(score(:,1));
pc1_cort = zscore(coeff(1:250,1));
pc1_gene = zscore(coeff(251:end,1));

save('./data/processed/main_decomp.mat','coeff','score','explained','GeneNames','seed_ind','used_seed_voxel_coords','used_seed_mni_coords','pc1_thal','pc1_cort','pc1_gene')

for i = 1:3
writematrix(score(:,i),['./data/processed/PC',num2str(i),'_thal.txt'],'Delimiter',' ')
end

SeedDists = squareform(pdist(seed_coords(logical(seed_ind),:)));
writematrix(SeedDists,'SeedDists.txt','Delimiter',' ')