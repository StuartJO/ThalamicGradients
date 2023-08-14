
TractGeneData = load('CompiledTractGeneData.mat');
CompiledTractData = PrepareTractGeneData(TractGeneData);

save('./data/processed/TractGeneNorm.mat','-struct','CompiledTractData')

TractData_GeneData_norm = [CompiledTractData.TractData_norm CompiledTractData.GeneData_norm];

Ncorts = size(CompiledTractData.TractData_norm,2);

decomp = RunPCADecomp(TractData_GeneData_norm,CompiledTractData.seed_ind,1:Ncorts,Ncorts+1:size(TractData_GeneData_norm,2));

save('./data/processed/main_decomp.mat','-struct','decomp')

for i = 1:3
writematrix(decomp.score(:,i),['./data/processed/PC',num2str(i),'_thal.txt'],'Delimiter',' ')
end

SeedDists = squareform(pdist(seed_vox_coords(logical(seed_ind),:)));
writematrix(SeedDists,'SeedDists.txt','Delimiter',' ')

%%
TractGeneData = load('CompiledTractGeneData_AllGeneSeed.mat');
CompiledTractData = PrepareTractGeneData(TractGeneData);
save('./data/processed/TractGeneNorm_AllGeneSeed.mat','-struct','CompiledTractData')
TractData_GeneData_norm = [CompiledTractData.TractData_norm CompiledTractData.GeneData_norm];
Ncorts = size(CompiledTractData.TractData_norm,2);
decomp = RunPCADecomp(TractData_GeneData_norm,CompiledTractData.seed_ind,1:Ncorts,Ncorts+1:size(TractData_GeneData_norm,2));
save('./data/processed/decomp_AllGeneSeed.mat','-struct','decomp')

TractGeneData = load('CompiledTractGeneData_Scha400.mat');
CompiledTractData = PrepareTractGeneData(TractGeneData);
save('./data/processed/TractGeneNorm_Scha400.mat','-struct','CompiledTractData')
TractData_GeneData_norm = [CompiledTractData.TractData_norm CompiledTractData.GeneData_norm];
Ncorts = size(CompiledTractData.TractData_norm,2);
decomp = RunPCADecomp(TractData_GeneData_norm,CompiledTractData.seed_ind,1:Ncorts,Ncorts+1:size(TractData_GeneData_norm,2));
save('./data/processed/decomp_Scha400.mat','-struct','decomp')

%% Run mouse PCA

RunAMBAdecomp

%% Run alternative decompositions

[AltEmbedding,DecompName] = RunAltDecomps();

save('./data/processed/alt_decomps','AltEmbedding','DecompName')
