function [AltEmbedding,DecompName] = RunAltDecomps()

%% Load in data

load('./data/preprocessed/CompiledTractGeneData.mat')

% Load in gene names and IDs for the AHBA data. These lists were obtained
% from Burt et al., 2020 Nature Neuroscience
EntrezIDs_all = dlmread('./data/preprocessed/AHBAEntrez.txt');
GeneIDs_burt_table = readtable('./data/preprocessed/AHBAgeneSymbol.txt','ReadVariableNames',false);

GeneIDs_burt = GeneIDs_burt_table.Var1;

% Find the names of the genes that existed in Gryglewski et al and Burt et
% al
[~,GeneIDs_burt_Gryglewski_ind] = ismember(GenesKept,EntrezIDs_all);

GeneIDs_valid = GeneIDs_burt(GeneIDs_burt_Gryglewski_ind);

% Load in the Phillips et al genes

PhillipsGeneTable = readtable('./data/preprocessed/PhillipsMouseThalGenes.xlsx');

PhillipsMouseGenes = PhillipsGeneTable.GeneSymbol;

PhillipsMouseGenesUpper = upper(PhillipsMouseGenes);

PhillipsMouse_HumanGenesOverlap_idx = find(ismember(GeneIDs_valid, PhillipsMouseGenesUpper));

% Normalise tract and gene data
TractData_norm = BF_NormalizeMatrix(ThalSeedAvg(:,1:250),'scaledSigmoid');
GeneData_norm = BF_NormalizeMatrix(SeedGene_kept,'scaledSigmoid');

TractData_GeneData_norm = [TractData_norm GeneData_norm];

%% Run joint decomposition with diffusion embedding
       
data = TractData_GeneData_norm';
sparse_data = data;
sparse_data(data < prctile(data,90)) = 0; 
cosine_similarity = 1-squareform(pdist(sparse_data','cosine'));

if ~all(conncomp(graph(abs(cosine_similarity),'lower')) == 1) 
    error('Graph is not connected.')
end

[AltEmbedding{1}] = diffusion_mapping(cosine_similarity, length(cosine_similarity), 0.5, 0);

% Just so happens with our data the direction of the gradient is reversed.
% because the direction is arbitary, we can reverse it by multiplying it by
% -1
AltEmbedding{1} = AltEmbedding{1}.*-1;

%% Run joint decomposition with Phillips mouse genes with PCA
Tract_PhillipsGenes = [TractData_norm GeneData_norm(:,PhillipsMouse_HumanGenesOverlap_idx)];

[~,AltEmbedding{2}] = pca(Tract_PhillipsGenes);

%% Run decomposition with tract and gene data seperately using PCA

[~,AltEmbedding{3}] = pca(TractData_norm);

[~,AltEmbedding{4}] = pca(GeneData_norm);

%% Run the decomposition of the first 10 PCs from the seperate tract and gene decompositions

[~,AltEmbedding{5}] = pca(zscore([AltEmbedding{3}(:,1:10) AltEmbedding{4}(:,1:10)]));

%% Averaged cosine similarity, DE 

simmat{1}= TractData_norm;
simmat{2} = GeneData_norm;

affinity_matrix = cell(1,2);

for i = 1:2
data = simmat{i}';
sparse_data = data;
sparse_data(data < prctile(data,90)) = 0; 
affinity_matrix{i} = 1-squareform(pdist(sparse_data','cosine'));
end
meanAffMat = (affinity_matrix{1}+affinity_matrix{2})./2;

if ~all(conncomp(graph(abs(meanAffMat),'lower')) == 1) 
    error('Graph is not connected.')
end

[AltEmbedding{6}] = diffusion_mapping(meanAffMat, length(affinity_matrix), 0.5, 0);

%% 

DecompName = {'Diffusion embedding','PCA on connectivity and Phillips genes','PCA on connectivity data only','PCA on gene-expression data only','PCA on first 10 connectvity and gene PCs','Diffusion embedding on mean affinity matrix'};
