
%% RunAnalysis_1

disp('Performing all decompositions')

addpath(genpath('./'))

%GetNiftiData()

TractGeneData = load('./data/preprocessed/CompiledTractGeneData_rand500.mat');
CompiledTractData = PrepareTractGeneData(TractGeneData);

save('./data/processed/TractGeneNorm_rand500.mat','-struct','CompiledTractData')

TractData_GeneData_norm = [CompiledTractData.TractData_norm CompiledTractData.GeneData_norm];

Ncorts = size(CompiledTractData.TractData_norm,2);

decomp = RunPCADecomp(TractData_GeneData_norm,CompiledTractData.seed_ind,1:Ncorts,Ncorts+1:size(TractData_GeneData_norm,2));
decomp.GeneNames_human = CompiledTractData.GeneNames_human;
save('./data/processed/decomp_rand500.mat','-struct','decomp')

for i = 1:3
writematrix(decomp.score(:,i),['./data/processed/PC',num2str(i),'_thal.txt'],'Delimiter',' ')
end

load('./data/preprocessed/MNI_Seed_voxelData.mat','seed_vox_coords')

SeedDists = squareform(pdist(seed_vox_coords(logical(TractGeneData.seed_ind),:)));
writematrix(SeedDists,'./data/processed/SeedDists.txt','Delimiter',' ')

%%
TractGeneData = load('./data/preprocessed/CompiledTractGeneData_AllGeneSeed.mat');
CompiledTractData = PrepareTractGeneData(TractGeneData);
save('./data/processed/TractGeneNorm_AllGeneSeed.mat','-struct','CompiledTractData')
TractData_GeneData_norm = [CompiledTractData.TractData_norm CompiledTractData.GeneData_norm];
Ncorts = size(CompiledTractData.TractData_norm,2);
decomp = RunPCADecomp(TractData_GeneData_norm,CompiledTractData.seed_ind,1:Ncorts,Ncorts+1:size(TractData_GeneData_norm,2));
decomp.GeneNames_human = CompiledTractData.GeneNames_human;
save('./data/processed/decomp_AllGeneSeed.mat','-struct','decomp')

TractGeneData = load('./data/preprocessed/CompiledTractGeneData_Scha400.mat');
CompiledTractData = PrepareTractGeneData(TractGeneData);
save('./data/processed/TractGeneNorm_Scha400.mat','-struct','CompiledTractData')
TractData_GeneData_norm = [CompiledTractData.TractData_norm CompiledTractData.GeneData_norm];
Ncorts = size(CompiledTractData.TractData_norm,2);
decomp = RunPCADecomp(TractData_GeneData_norm,CompiledTractData.seed_ind,1:Ncorts,Ncorts+1:size(TractData_GeneData_norm,2));
decomp.GeneNames_human = CompiledTractData.GeneNames_human;
save('./data/processed/decomp_Scha400.mat','-struct','decomp')

%% Run mouse PCA

RunAMBAdecomp

%% Run alternative decompositions

[AltEmbedding,DecompName] = RunAltDecomps();

save('./data/processed/alt_decomps','AltEmbedding','DecompName')