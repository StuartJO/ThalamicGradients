function RunMainDecomp()

load('./data/processed/TractGeneNorm.mat')

TractData_GeneData_norm = [TractData_norm GeneData_norm];

decomp = struct;

decomp.input = TractData_GeneData_norm;

[decomp.coeff,decomp.score,~,~,decomp.explained] = pca(TractData_GeneData_norm);

load('./data/ancillary/MNI_Seed_voxelData.mat','seed_vox_coords','seed_mni_coords')

decomp.used_seed_voxel_coords = used_seed_voxel_coords;
decomp.used_seed_mni_coords = used_seed_mni_coords;

decomp.pcs_thal = zscore(decomp.score(:,:));
decomp.pcs_cort = zscore(decomp.coeff(1:250,:));
decomp.pcs_gene = zscore(decomp.coeff(251:end,:));

save('./data/processed/main_decomp.mat','-struct','decomp')

for i = 1:3
writematrix(decomp.score(:,i),['./data/processed/PC',num2str(i),'_thal.txt'],'Delimiter',' ')
end

SeedDists = squareform(pdist(seed_vox_coords(logical(seed_ind),:)));
writematrix(SeedDists,'SeedDists.txt','Delimiter',' ')