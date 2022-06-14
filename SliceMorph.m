figure('Position',[0 0 1920 1080])

slice_xlimits = [98.0211  108.5175];
slice_ylimits = [93.4533  127.2028];
slice_zlimits = [116.0609  143.8561];
load('MNI_Seed_voxelData.mat')
brain = niidata;
load('MNI_1mm_mask2.mat','brain_mask')
z_loc = 82;
brain_rescaled = rescale(brain,1,256);
brain_rescaled(~brain_mask)=NaN;
brain_slice = brain_rescaled(:,:,z_loc)';
colormap(gray(256))
slicesurf = surf(repmat(z_loc, size(brain_slice)), brain_slice);
slicesurf.EdgeAlpha=0;
%slicesurf.FaceAlpha=0;
slicesurf.FaceLighting='none';
slicesurf.Clipping='off';

hold on
thalCmap=turbo(256);
[newslicevals,newCmap] = ProjThalData2Slice(pc_thal,z_loc,seed_voxel_coords,thalCmap,1.75);
% 
z_loc = 82;
brain_slice = brain_rescaled(:,:,z_loc)';
colormap(gray(256))
slicesurf2 = surf(repmat(z_loc+.1, size(brain_slice)), newslicevals);
slicesurf2.EdgeAlpha=0;
slicesurf2.FaceAlpha=0;
slicesurf2.FaceLighting='none';
slicesurf2.Clipping='off';
slicesurf2.CData=newslicevals;
colormap(newCmap)


thalCmap=turbo(256);
maskslice = newslicevals;
% 
maskslice(maskslice>256) = NaN;
z_loc = 82;
brain_slice = brain_rescaled(:,:,z_loc)';
sliceMASK = surf(repmat(87, size(brain_slice)), maskslice);
sliceMASK.EdgeAlpha=0;
sliceMASK.FaceAlpha=0;
sliceMASK.FaceLighting='none';
sliceMASK.Clipping='off';
sliceMASK.CData=maskslice;

view([0 90])

axis off
axis tight
axis equal
axis vis3d

xlim(slice_xlimits)
ylim(slice_ylimits)
zlim(slice_zlimits)

pc_thal = (zscore(score{3}(:,1)'));
pc_patch_coords_start = linspace(121.75,92,922)';
cmap = turbo(256);
pc_thal_colors_ind = ceil(rescale(pc_thal,1,size(cmap,1)));
pc_thal_colors_ind(isnan(pc_thal)) = 1;
pc_thal_colors = cmap(pc_thal_colors_ind,:);


x_coords = [121.75+5 121.75+10 121.75+10 121.75+5];
z_coords = [85 85 85 85];


seed_voxel_coords = seeds_vox(logical(seed_ind),:);
dist2slice = pdist2(seed_voxel_coords(:,3),z_loc);
minDist = min(dist2slice);
seeds2plot_ind = dist2slice==minDist;

vox_spacing = .2;
r = linspace(0,1,30);
zoffset = 2;

seeds2plot = seed_voxel_coords(seeds2plot_ind,:);

seeds2plot_colors = pc_thal_colors(seeds2plot_ind,:);
[x,y,z]=sphere(64);
sphere_mesh = surf2patch(x,y,z,'triangles');
for i = 1:length(seeds2plot)
    sphere_mesh_roi = sphere_mesh;
    sphere_mesh_roi.vertices = (sphere_mesh.vertices.*.25) + [seeds2plot(i,1:2)+.5 82.5];
    sSeeds(i) = patch(sphere_mesh_roi,'EdgeColor','none','FaceColor',[0 0 0],'Clipping','off','FaceLighting','none');
    material dull     
end

print(['./NewGIFs/SliceMorph_START.png'],'-dpng')

for i = 1:921
y_coords = [pc_patch_coords_start(i) pc_patch_coords_start(i) pc_patch_coords_start(i+1) pc_patch_coords_start(i+1)]+5; 
f = [1 2 3 4];
v = [x_coords(1) y_coords(1) z_coords(1);x_coords(2) y_coords(2) z_coords(2);x_coords(3) y_coords(3) z_coords(3);x_coords(4) y_coords(4) z_coords(4)];

pc_patches(i) = patch('Faces',f,'Vertices',v,'FaceColor',pc_thal_colors(i,:),'Clipping','off','EdgeColor','none','FaceLighting','none');
end


for j = 1:30
    for i = 1:921
    y_coords = [pc_patch_coords_start(i) pc_patch_coords_start(i) pc_patch_coords_start(i+1) pc_patch_coords_start(i+1)]+5; 
    v = [x_coords(1) y_coords(1) z_coords(1);x_coords(2) y_coords(2) z_coords(2);x_coords(3) y_coords(3) z_coords(3);x_coords(4) y_coords(4) z_coords(4)];

    voxel_coord = seed_voxel_coords(i,:)+[.5 .5 0];
    
    v2 = [voxel_coord+[-vox_spacing -vox_spacing zoffset]; voxel_coord+[vox_spacing -vox_spacing zoffset];voxel_coord+[vox_spacing vox_spacing zoffset];voxel_coord+[-vox_spacing vox_spacing zoffset]];
    
    v3 = find_point_on_line(v,v2,r(j));
    
    pc_patches(i).Vertices = v3;
    
    if seeds2plot_ind(i)==0
        pc_patches(i).FaceAlpha=1-r(j);
    end
    
    end
%     disp(num2str(j))
     pause(.1)
    print(['./NewGIFs/SliceMorph1_',num2str(j),'.png'],'-dpng')
end

sliceMASK.FaceAlpha=0;

% seeds2plot = seed_voxel_coords(seeds2plot_ind,:);
% 
% seeds2plot_colors = pc_thal_colors(seeds2plot_ind,:);
% 
% for i = 1:length(seeds2plot)
%     sSeeds(i).FaceColor = seeds2plot_colors(i,:);
%     material dull     
% end

%vox_spacing2 = .875;
vox_spacing2 = 1.75;
r2 = linspace(0,1,30);
r3 = [zeros(1,20) linspace(0,1,10)];
    for i = 1:length(seeds2plot)
        sSeeds(i).FaceAlpha = 1;   
    end
for j = 1:30
    for i = 1:921
    voxel_coord = seed_voxel_coords(i,:)+[.5 .5 0];
    
    v = [voxel_coord+[-vox_spacing -vox_spacing zoffset]; voxel_coord+[vox_spacing -vox_spacing zoffset];voxel_coord+[vox_spacing vox_spacing zoffset];voxel_coord+[-vox_spacing vox_spacing zoffset]];
    
    v2 = [voxel_coord+[-vox_spacing2 -vox_spacing2 zoffset]; voxel_coord+[vox_spacing2 -vox_spacing2 zoffset];voxel_coord+[vox_spacing2 vox_spacing2 zoffset];voxel_coord+[-vox_spacing2 vox_spacing2 zoffset]];
        
    v3 = find_point_on_line(v,v2,r2(j));
    
    pc_patches(i).Vertices = v3;
    
    if seeds2plot_ind(i)==0
        pc_patches(i).FaceAlpha=0;
    else
        pc_patches(i).FaceAlpha=1-r3(j);
    end
    
    end
    
    if j == 10
    for i = 1:length(seeds2plot)
        sSeeds(i).FaceAlpha = 0;   
    end
    end
    
    slicesurf2.FaceAlpha = r3(j);
%     disp(num2str(j))
     pause(.1)
     print(['./NewGIFs/SliceMorph2_',num2str(j),'.png'],'-dpng')
end
