function decomp = RunPCADecomp(Input, seed_ind, TractInds, GeneInds)
% Run PCA decomposition on input data and return decomposition results.
%
% Inputs:
%   Input - Input data matrix (rows: samples, columns: features)
%   seed_ind - Logical vector indicating whether each sample corresponds to a seed
%   TractInds - Indices of columns corresponding to tract-related features
%   GeneInds - Indices of columns corresponding to gene-related features
%
% Outputs:
%   decomp - Structure containing PCA decomposition results and related information
%     .input - Original input data matrix
%     .coeff - Coefficients of the principal components
%     .score - PCA scores (projection of data onto principal components)
%     .explained - Variance explained by each principal component
%     .used_seed_voxel_coords - Voxel coordinates of used seed samples
%     .used_seed_mni_coords - MNI coordinates of used seed samples
%     .pcs_thal - Standardized PCA scores for thalamic data
%     .pcs_cort - Standardized PCA scores for cortical tract data
%     .pcs_gene - Standardized PCA scores for gene data

decomp = struct;

% Store the original input data
decomp.pca_input = Input;

% Perform PCA decomposition on the input data
[decomp.coeff, decomp.score, ~, ~, decomp.explained] = pca(Input);

% Load seed voxel coordinates and MNI coordinates
load('./data/preprocessed/MNI_Seed_voxelData.mat', 'seed_vox_coords', 'seed_mni_coords')

% Store the used seed voxel coordinates and MNI coordinates
decomp.used_seed_voxel_coords = seed_vox_coords(logical(seed_ind), :);
decomp.used_seed_mni_coords = seed_mni_coords(logical(seed_ind), :);

% Standardize and store the PCA scores for thalamic data
decomp.pcs_thal = zscore(decomp.score(:, :));

% Standardize and store the PCA scores for cortical tract data
decomp.pcs_cort = zscore(decomp.coeff(TractInds, :));

% Standardize and store the PCA scores for gene data
decomp.pcs_gene = zscore(decomp.coeff(GeneInds, :));
