function GetCorticalSpinTestPerms()

Niters = 10000;

lh_sphere = read_surface('lh.sphere');
load('fsaverage_surface_data.mat','lh_rand500')

centroid = zeros(250,3);

for i = 1:250
    centroid(i,:)=mean(lh_sphere.vertices(lh_rand500==i,:));
end

perm_id = rotate_parcellation_single(centroid,Niters);
save('./data/preprocessed/CorticalSpinTestPerms.mat','perm_id')