function RunAMBAdecomp()

% Load in the mouse connectivity data
% As per the description here:
% https://github.com/benfulcher/mouseStructuralConnectivity/blob/master/ImportConnectomeData.m,
% Conn_W{1} is the connectivity for within a single (right) hemisphere.
% While we are at it lets also get the names of each region (RegionStruct)

load ('./data/preprocessed/Mouse_Connectivity_Data.mat','Conn_W','RegionStruct')

% Load in the AMBA gene data

load ('./data/preprocessed/AllenGeneDataset_19419.mat','GeneExpData','geneInfo','structInfo')

GeneData_fullmouse = GeneExpData.comb.energy;

% Get the list of Phillips genes
PhillipsGeneTable = readtable('PhillipsMouseThalGenes.xlsx');

GeneIDs_phillips = PhillipsGeneTable.GeneSymbol;

ThalRegions = 88:122;
%ThalRegions = [88:93 95:99 101:117 121 122];

GeneData_thal = GeneData_fullmouse(ThalRegions,:);

% Get only genes for which we have complete data

NanGenes = isnan(sum(GeneData_thal,1));

GeneData_thal_complete = GeneData_thal(:,~NanGenes);

GeneIDs_thal_complete = geneInfo.acronym(~NanGenes);

% Get only genes which overlap with Phillips et al

GeneIDs_thal_complete_phillips_idx = find(ismember(upper(GeneIDs_thal_complete), upper(GeneIDs_phillips)));

mouse_GeneNames = GeneIDs_thal_complete(GeneIDs_thal_complete_phillips_idx);

% Extract connectivity and expression for thalamic regions

CortRegions = 1:38;

mouse_GeneData_norm = BF_NormalizeMatrix(GeneData_thal_complete(:,GeneIDs_thal_complete_phillips_idx),'scaledSigmoid');

mouse_TractData = Conn_W{1};

mouse_TractData_norm = BF_NormalizeMatrix(mouse_TractData(ThalRegions,CortRegions),'scaledSigmoid');

mouse_TractData_GeneData_norm = [mouse_TractData_norm mouse_GeneData_norm];

[mouse_coeff,mouse_score,~,~,mouse_explained] = pca(mouse_TractData_GeneData_norm);

mouse_pcs_thal = zscore(mouse_score);
mouse_pcs_cort = zscore(mouse_coeff(1:38,:));
mouse_pcs_gene = zscore(mouse_coeff(39:end,:));

save('./data/processed/mouse_decomp.mat','mouse_coeff','mouse_score','mouse_explained','mouse_GeneNames','mouse_pcs_thal','mouse_pcs_cort','mouse_pcs_gene')