function [mouse_CortFeatures,mouse_ThalHierarchy] = GetAMBAcortexdata()
   
load ('./data/preprocessed/AllenGeneDataset_19419.mat','structInfo')

AMBA_cort_data = load('./data/preprocessed/Data_AMBAcortex.mat');

MouseHierarchy = readtable('./data/preprocessed/TCCT_CCconf_iter.xls','Range','A2:E62');

[~,idx] = ismember(lower(structInfo.acronym(1:38)),lower(MouseHierarchy.Var2));

mouse_CortHierarchy = nan(38,1);

NotNan = find(idx);

mouse_CortHierarchy(NotNan) = MouseHierarchy.Var5(idx(NotNan));

[~,~,iy] = intersect(lower(structInfo.acronym),lower(AMBA_cort_data.structInfo.acronym),'stable');

mouse_CortFeatures.data = AMBA_cort_data.dataMatrix(iy,:);

mouse_CortFeatures.data(:,6) = mouse_CortHierarchy;

mouse_CortFeatures.name = {'Grik2 expression','Pvalb expression','Grin3a expression','PV cell density','Axonal input strength','Hierarchical level','Cytoarchitecture type','T1w:T2w','PC1 transcription'};

ThalRegions = 88:122;

[~,idx] = ismember(lower(structInfo.acronym(ThalRegions)),lower(MouseHierarchy.Var2));

mouse_ThalHierarchy = nan(length(ThalRegions),1);

NotNan = find(idx);

mouse_ThalHierarchy(NotNan) = MouseHierarchy.Var5(idx(NotNan));
