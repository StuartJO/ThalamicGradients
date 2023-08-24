function [AltEmbedding, DecompName] = RunAltDecomps()
% This function runs alternative decompositions on tract and gene data

%% Load in data

% Load compiled tract and gene data
load('./data/preprocessed/CompiledTractGeneData_Rand500.mat','GenesKept','ThalSeedAvg','ThalSeedGenesKept')

% Load gene names and IDs from Burt et al., 2020 Nature Neuroscience
EntrezIDs_all = dlmread('./data/preprocessed/AHBAEntrez.txt');
GeneIDs_burt_table = readtable('./data/preprocessed/AHBAgeneSymbol.txt', 'ReadVariableNames', false);
GeneIDs_burt = GeneIDs_burt_table.Var1;

% Find gene names existing in both Gryglewski et al. and Burt et al.
[~, GeneIDs_burt_Gryglewski_ind] = ismember(GenesKept, EntrezIDs_all);
GeneIDs_valid = GeneIDs_burt(GeneIDs_burt_Gryglewski_ind);

% Load Phillips et al. mouse genes
PhillipsGeneTable = readtable('./data/preprocessed/PhillipsMouseThalGenes.xlsx');
PhillipsMouseGenes = PhillipsGeneTable.GeneSymbol;
PhillipsMouseGenesUpper = upper(PhillipsMouseGenes);
PhillipsMouse_HumanGenesOverlap_idx = ismember(GeneIDs_valid, PhillipsMouseGenesUpper);

% Normalize tract and gene data
TractData_norm = BF_NormalizeMatrix(ThalSeedAvg(:, 1:250), 'scaledSigmoid');
GeneData_norm = BF_NormalizeMatrix(ThalSeedGenesKept, 'scaledSigmoid');
TractData_GeneData_norm = [TractData_norm GeneData_norm];

AltEmbedding = cell(1,6);

%% Run joint decomposition with diffusion embedding

data = TractData_GeneData_norm';
sparse_data = data;
sparse_data(data < prctile(data, 90)) = 0;
cosine_similarity = 1 - squareform(pdist(sparse_data', 'cosine'));

if ~all(conncomp(graph(abs(cosine_similarity), 'lower')) == 1)
    error('Graph is not connected.')
end

[AltEmbedding{1}] = diffusion_mapping(cosine_similarity, length(cosine_similarity), 0.5, 0);
AltEmbedding{1} = AltEmbedding{1} .* -1; % Reverse the gradient direction

%% Run joint decomposition with Phillips mouse genes using PCA

Tract_PhillipsGenes = [TractData_norm GeneData_norm(:, PhillipsMouse_HumanGenesOverlap_idx)];
[~, AltEmbedding{4}] = pca(Tract_PhillipsGenes);

%% Run separate decomposition with PCA for tract and gene data

[~, AltEmbedding{3}] = pca(TractData_norm);
[~, AltEmbedding{2}] = pca(GeneData_norm);

%% Run decomposition of the first 10 PCs from separate tract and gene decompositions

[~, AltEmbedding{5}] = pca(zscore([AltEmbedding{2}(:, 1:10) AltEmbedding{3}(:, 1:10)]));

%% Calculate averaged cosine similarity and diffusion embedding

simmat{1} = TractData_norm;
simmat{2} = GeneData_norm;

affinity_matrix = cell(1, 2);

for i = 1:2
    data = simmat{i}';
    sparse_data = data;
    sparse_data(data < prctile(data, 90)) = 0;
    affinity_matrix{i} = 1 - squareform(pdist(sparse_data', 'cosine'));
end
meanAffMat = (affinity_matrix{1} + affinity_matrix{2}) ./ 2;

if ~all(conncomp(graph(abs(meanAffMat), 'lower')) == 1)
    error('Graph is not connected.')
end

[AltEmbedding{6}] = diffusion_mapping(meanAffMat, length(affinity_matrix), 0.5, 0);

%% Assign decomposition names

DecompName = {'Diffusion embedding', 'PCA on gene-expression data only', 'PCA on connectivity data only', 'PCA on connectivity and Phillips genes', 'PCA on first 10 connectivity and gene PCs', 'Diffusion embedding on mean affinity matrix'};
