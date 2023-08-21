function GetNiftiData()

% This function extracts relevant nifti data and gets the seed coordinates 
% and puts them into a .mat file

nifti_info = niftiinfo('brain_1mm.nii.gz');
braindata = niftiread('brain_1mm.nii.gz');

mm = '1.75';

seed_mni_coords = dlmread([DIR,'/seeds_',mm,'mm_ind.txt']);