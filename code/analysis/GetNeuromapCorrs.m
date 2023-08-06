function neuromap_corrs = GetNeuromapCorrs(pc_cort)

if nargin < 1
    load('./data/processed/main_decomp.mat','pc1_cort')
    pc_cort = pc1_cort;
end

load('fsaverage_surface_data.mat','lh_rand500')

parc_size = 250;

files = dir('./data/preprocessed/fsaverage164k_annots/*.gii');
Nmaps = length(files);

neuromap = zeros(length(lh_rand500),Nmaps);

for i = 1:Nmaps
    g = gifti(files(i).name);
    neuromap(:,i) = g.cdata;
end

neuromap_corrs.neuromap_parc = zeros(parc_size,Nmaps);

for i = 1:parc_size
    neuromap_corrs.neuromap_parc(i,:) = nanmean(neuromap(lh_rand500==i,:),1);
end

load('./data/processed/CorticalSpinTestPerms.mat','perm_id')

neuromap_corrs.data = pc_cort;

neuromap_corrs.p_perm = zeros(Nmaps,1);

for i = 1:Nmaps
    neuromap_corrs.p_perm(i) = perm_sphere_p(neuromap_corrs.data,neuromap_corrs.neuromap_parc(:,i),perm_id,'Pearson');
end

neuromap_corrs.corr = corr(neuromap_corrs.data,neuromap_corrs.neuromap_parc,'Rows','complete')';

neuromap_corrs.corr_sig = neuromap_corrs.p_perm<.05;

NeuroMaps_names = readtable('NeuroMaps_names.xlsx');
for i = 1:size(NeuroMaps_names,1)
    Index = find(contains(NeuroMaps_names.Filename,files(i).name));
    neuromap_corrs.name{i} = NeuroMaps_names.ShortName{Index};
    neuromap_corrs.description{i} = NeuroMaps_names.Description{Index};
end