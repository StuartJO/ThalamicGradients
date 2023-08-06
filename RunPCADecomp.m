function decomp = RunPCADecomp(Input,seed_ind,TractInds,GeneInds)

decomp = struct;

decomp.input = Input;

[decomp.coeff,decomp.score,~,~,decomp.explained] = pca(Input);

load('./data/ancillary/MNI_Seed_voxelData.mat','seed_vox_coords','seed_mni_coords')

decomp.used_seed_voxel_coords = seed_vox_coords(logical(seed_ind),:);
decomp.used_seed_mni_coords = seed_mni_coords(logical(seed_ind),:);

decomp.pcs_thal = zscore(decomp.score(:,:));
decomp.pcs_cort = zscore(decomp.coeff(TractInds,:));
decomp.pcs_gene = zscore(decomp.coeff(GeneInds,:));