function GetNiftiData()

% This function extracts relevant nifti data and gets the seed coordinates 
% and puts them into a .mat file

nifti_info = niftiinfo('./data/preprocessed/brain_1mm.nii.gz');
braindata = niftiread('./data/preprocessed/brain_1mm.nii.gz');

thalmask = niftiread('./data/preprocessed/Left_Thal.nii.gz');

mm = '1.75';

seeds = readtable(['./data/preprocessed/seeds_',mm,'mm.txt']);

seed_mni_coords = [seeds.Pos_x seeds.Pos_y seeds.Pos_z];

seed_vox_coords = ConvertNiftiCoords(nifti_info, seed_mni_coords, 'mm');

save('MNI_Seed_voxelData.mat','seed_vox_coords','seed_mni_coords','thalmask','braindata')