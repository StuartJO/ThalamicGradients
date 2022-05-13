
pc_patch_coords_start = linspace(130,40,922)';
cmap = turbo(256);
pc_thal_colors_ind = ceil(rescale(pc_thal,1,size(cmap,1)));
pc_thal_colors_ind(isnan(pc_thal)) = 1;
pc_thal_colors = cmap(pc_thal_colors_ind,:);

x = [-200 -200 -200 -200];
y = [230 240 240 230];
for i = 1:921
z = [pc_patch_coords_start(i) pc_patch_coords_start(i) pc_patch_coords_start(i+1) pc_patch_coords_start(i+1)]; 
f = [1 2 3 4];
v = [x(1) y(1) z(1);x(2) y(2) z(2);x(3) y(3) z(3);x(4) y(4) z(4)];

pc_patches(i) = patch('Faces',f,'Vertices',v,'FaceColor',pc_thal_colors(i,:),'Clipping','off','EdgeColor','none');
end



z_loc = 82;
brain_slice = brain(:,:,z_loc)';
colormap(gray(256))
slicesurf = surf(repmat(z_loc, size(brain_slice)), brain_slice);
slicesurf.EdgeAlpha=0;
slicesurf.Clipping='off';
colormap(gray(256))
xlim(xlimits)
ylim(ylimits)
zlim(zlimits)
axis tight

repmat(z_loc, size(brain_slice));