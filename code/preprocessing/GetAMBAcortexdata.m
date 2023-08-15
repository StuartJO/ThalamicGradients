function [mouse_CortFeatures, mouse_ThalHierarchy, mouse_corticalRegions, mouse_thalamicRegion] = GetAMBAcortexdata()
% GETAMBACKACORTEXDATA Retrieves cortex and thalamus data from the Allen Mouse Brain Atlas
%
%   [mouse_CortFeatures, mouse_ThalHierarchy, mouse_corticalRegions, mouse_thalamicRegion] = GetAMBAcortexdata()
%
%   Outputs:
%   - mouse_CortFeatures: A structure containing cortex-related features extracted
%     from the Allen Mouse Brain Atlas data. It includes a numeric data matrix
%     and feature names.
%       - 'data': A numeric matrix where each row represents a specific cortical
%         region, and each column corresponds to a different feature. The number
%         of rows corresponds to the number of matched cortical regions from the
%         AMBA data, and the number of columns is equal to the number of defined
%         features (including hierarchy level and other attributes). If a feature
%         is not measured for a cortical region, it will have a value of 'NaN'.
%       - 'name': A cell array of strings providing the names of the features
%         corresponding to each column of the 'data' matrix.
%
%   - mouse_ThalHierarchy: A numeric array containing hierarchical level
%     information for thalamic regions of interest extracted from the Allen Mouse
%     Brain Atlas data. It provides hierarchy levels for the 35 thalamic regions.
%     If no data is available for that region, it has a value of 'NaN'.
%
%   - mouse_corticalRegions: A cell array of strings representing the acronyms
%     of the matched cortical regions in the Allen Mouse Brain Atlas.
%
%   - mouse_thalamicRegion: A cell array of strings representing the acronyms
%     of the thalamic regions of interest extracted from the Allen Mouse Brain Atlas.
%
%   This function retrieves cortex-related features and hierarchy level
%   information from the Allen Mouse Brain Atlas data for further analysis or
%   visualization. The outputs can be used to explore relationships between the
%   defined features and hierarchy levels within the cortex and thalamic regions.

% Load structural information about brain regions
load('./data/preprocessed/AllenGeneDataset_19419.mat', 'structInfo')

% Load AMBA cortex data
AMBA_cort_data = load('./data/preprocessed/Data_AMBAcortex.mat');

% Load mouse hierarchy data
MouseHierarchy = readtable('./data/preprocessed/TCCT_CCconf_iter.xls', 'Range', 'A2:E62');

% Store acronyms of matched cortical regions
mouse_corticalRegions = structInfo.acronym(1:38);

% Match cortex data to mouse hierarchy
[~, idx] = ismember(lower(structInfo.acronym(1:38)), lower(MouseHierarchy.Var2));
mouse_CortHierarchy = nan(38, 1);
NotNan = find(idx);
mouse_CortHierarchy(NotNan) = MouseHierarchy.Var5(idx(NotNan));

% Match cortex data to AMBA cortex data
[~, ~, iy] = intersect(lower(structInfo.acronym), lower(AMBA_cort_data.structInfo.acronym), 'stable');
mouse_CortFeatures.data = AMBA_cort_data.dataMatrix(iy, :);
mouse_CortFeatures.data(:, 6) = mouse_CortHierarchy;

% Define cortex feature names
mouse_CortFeatures.name = {'Grik2 expression', 'Pvalb expression', 'Grin3a expression', 'PV cell density', 'Axonal input strength', 'Hierarchical level', 'Cytoarchitecture type', 'T1w:T2w', 'PC1 transcription'};

% Define thalamic regions of interest
ThalRegions = 88:122;

% Store acronyms of thalamic regions
mouse_thalamicRegion = structInfo.acronym(ThalRegions);

% Match thalamic data to mouse hierarchy
[~, idx] = ismember(lower(structInfo.acronym(ThalRegions)), lower(MouseHierarchy.Var2));
mouse_ThalHierarchy = nan(length(ThalRegions), 1);
NotNan = find(idx);
mouse_ThalHierarchy(NotNan) = MouseHierarchy.Var5(idx(NotNan));
