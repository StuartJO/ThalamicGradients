function neuromap_corrs = GetNeuromapCorrs(pc_cort,parc)

% This function calculates correlations between cortical principal 
% components (PCs) and preprocessed neuromaps. Note the neuromaps are on
% the fsaverage 164k surface.

% Inputs:
%   pc_cort: A matrix representing cortical principal components (PCs).
%            Each row corresponds to a different data point, and each column
%            corresponds to a different PC.
%   parc:    An array representing the parcellation of the fsaverage 164k surface.
%            Each element corresponds to the parcel index of a vertex on the surface.

% Outputs:
%   neuromap_corrs: A structure containing the following fields:
%       - neuromap_parc: A matrix where each row corresponds to a parcel and each
%                        column corresponds to a neuroanatomical map. This matrix
%                        contains mean values of neuroanatomical maps within each parcel.
%       - data: The input cortical principal components matrix.
%       - p_perm: An array of permutation test p-values for each neuroanatomical map.
%       - corr: A matrix of Pearson correlations between PC data and neuroanatomical maps.
%       - corr_sig: A logical array indicating significant correlations (p < 0.05).
%       - name: Cell array of short names for each neuroanatomical map.
%       - description: Cell array of descriptions for each neuroanatomical map.

% Check if input argument is provided, if not, load default pc_cort
if nargin < 1
    load('./data/processed/decomp_rand500.mat', 'pc1_cort')
    pc_cort = pc1_cort;
end

if nargin < 2
    % Load preprocessed data for fsaverage surface 
    load('fsaverage_surface_data.mat', 'lh_rand500')
    parc = lh_rand500;
end

% Get the number of parcels. This assumes sequential numbering of parcels
parc_size = max(parc);

% Load preprocessed neuroanatomical map files
files = dir('./data/preprocessed/fsaverage164k_annots/*.gii');
Nmaps = length(files);

% Initialize matrix to store neuroanatomical maps
neuromap = zeros(length(parc), Nmaps);

% Load data from each file into the neuromap matrix
for i = 1:Nmaps
    g = gifti(files(i).name);
    neuromap(:, i) = g.cdata;
end

% Initialize a structure to store results of correlation analysis
neuromap_corrs.neuromap_parc = zeros(parc_size, Nmaps);

% Calculate mean values for each parcellated region of the neuroanatomical maps
for i = 1:parc_size
    neuromap_corrs.neuromap_parc(i, :) = nanmean(neuromap(parc == i, :), 1);
end

% Load preprocessed permutation IDs for statistical analysis
load('./data/processed/CorticalSpinTestPerms.mat', 'perm_id')

% Store input data in the results structure
neuromap_corrs.data = pc_cort;

% Initialize an array to store permutation test p-values
neuromap_corrs.p_perm = zeros(Nmaps, 1);

% Perform permutation test for each neuroanatomical map and store p-values
for i = 1:Nmaps
    neuromap_corrs.p_perm(i) = perm_sphere_p(neuromap_corrs.data, neuromap_corrs.neuromap_parc(:, i), perm_id, 'Pearson');
end

% Calculate Pearson correlations between PC data and neuroanatomical maps
neuromap_corrs.corr = corr(neuromap_corrs.data, neuromap_corrs.neuromap_parc, 'Rows', 'complete')';

% Identify significant correlations based on a threshold of 0.05
neuromap_corrs.corr_sig = neuromap_corrs.p_perm < 0.05;

% Read NeuroMaps_names from a table
NeuroMaps_names = readtable('NeuroMaps_names.xlsx');

% Assign NeuroMaps_names and descriptions to the results structure
for i = 1:size(NeuroMaps_names, 1)
    Index = find(contains(NeuroMaps_names.Filename, files(i).name));
    neuromap_corrs.name{i} = NeuroMaps_names.ShortName{Index};
    neuromap_corrs.description{i} = NeuroMaps_names.Description{Index};
end