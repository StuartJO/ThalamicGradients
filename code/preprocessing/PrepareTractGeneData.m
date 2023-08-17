function CompiledTractData = PrepareTractGeneData(TractGeneData)

% CompiledTractData contains the following variables:
% GenesKept = Entrez IDs of the genes obtained from Gryglewski
% SUB = IDs of HCP subjects whose tractography data was used
% SeedGene_kept = matrix of seeds-by-gene-expression
% ThalSeed = matrix of seeds-by-ROI for each SUB.
% ThalSeedAvg = ThalSeed averaged (mean) across subjects
% seed_ind = index of which of the 1811 seeds were kept

% Load in gene names and IDs for the AHBA data. These lists were obtained
% from Burt et al., 2020 Nature Neuroscience
EntrezIDs_all = dlmread('AHBAEntrez.txt');
GeneIDs_burt_table = readtable('AHBAgeneSymbol.txt','ReadVariableNames',false);

GeneIDs_burt = GeneIDs_burt_table.Var1;

% Find the names of the genes that existed in Gryglewski et al and Burt et
% al (i.e., the overlap)
[~,GeneIDs_burt_Gryglewski_ind] = ismember(TractGeneData.GenesKept,EntrezIDs_all);

GeneNames = GeneIDs_burt(GeneIDs_burt_Gryglewski_ind);

TractData_norm = BF_NormalizeMatrix(TractGeneData.ThalSeedAvg,'scaledSigmoid');
GeneData_norm = BF_NormalizeMatrix(TractGeneData.SeedGene_kept,'scaledSigmoid');

load('./data/ancillary/MNI_Seed_voxelData.mat','seed_vox_coords','seed_mni_coords')

used_seed_voxel_coords = seed_vox_coords(logical(TractGeneData.seed_ind),:);
used_seed_mni_coords = seed_mni_coords(logical(TractGeneData.seed_ind),:);

CompiledTractData.TractData_norm = TractData_norm;
CompiledTractData.GeneData_norm = GeneData_norm;
CompiledTractData.seed_ind = TractGeneData.seed_ind;
CompiledTractData.used_seed_voxel_coords = used_seed_voxel_coords;
CompiledTractData.used_seed_mni_coords = used_seed_mni_coords;
CompiledTractData.GeneNames_human = GeneNames;