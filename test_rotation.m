
sSeeds = scatter3(seed_voxel_coords(:,1),seed_voxel_coords(:,2),seed_voxel_coords(:,3),60,'filled','k');

% sSeeds = scatter3(seed_voxel_coords(:,1),seed_voxel_coords(:,2),seed_voxel_coords(:,3),60,pc_thal,'filled');
sSeeds.Clipping='off';
% sSeeds.SizeData=60;
% colormap(turbo(256))

ax = gca;
origin = sum([get(ax,'xlim')' get(ax,'ylim')' get(ax,'zlim')'])/2;

pc_line_coords_start = [ones(921,1)*142 ones(921,1)*145 linspace(100,80,921)'];
pc_line_coords_end = [ones(921,1)*142 ones(921,1)*155 linspace(100,80,921)'];

pc_line_coords_start_rot = rotate_coords(pc_line_coords_start,views(end,:),1,origin);
pc_line_coords_end_rot = rotate_coords(pc_line_coords_end,views(end,:),1,origin);

cmap = turbo(256);

 pc_thal_colors_ind = ceil(rescale(pc_thal,1,size(cmap,1)));
 pc_thal_colors_ind(isnan(pc_thal)) = 1;
 pc_thal_colors = cmap(pc_thal_colors_ind,:);
for i = 1:921
    x = [pc_line_coords_start(i,1) pc_line_coords_end(i,1)];
    y = [pc_line_coords_start(i,2) pc_line_coords_end(i,2)];
    z = [pc_line_coords_start(i,3) pc_line_coords_end(i,3)];
    l_score(i) = plot3(x,y,z,'Color',pc_thal_colors(i,:)); 
    l_score(i).Clipping = 'off';
    hold on
end
view(views(end,:))
% ylim(ylimits)
% zlim(zlimits)
% axis tight
for i = 1:921
% direction = [0 0 1];
% rotate(l_score(i),direction,46)

direction = [0 -1 0];
rotate(l_score(i),direction,-13)
end

    view(views(end,:))
    xlim(xlimits_range(end,:));
    ylim(ylimits_range(end,:));
    zlim(zlimits_range(end,:));
    axis tight
axis equal
axis vis3d
   hold on 
pc_patch_coords_start = [ones(922,1)*142 ones(922,1)*145 linspace(100,80,922)'];
pc_patch_coords_end = [ones(922,1)*142 ones(922,1)*155 linspace(100,80,922)'];
for i = 1:921
x = [142 142 142 142];
y = [145 145 155 155];
z = [pc_patch_coords_start(i,3) pc_patch_coords_start(i+1,3) pc_patch_coords_start(i+1,3) pc_patch_coords_start(i,3)]; 
f = [1 2 3 4];
v = [x(1) y(1) z(1);x(2) y(2) z(2);x(3) y(3) z(3);x(4) y(4) z(4)];

pc_patches(i) = patch('Faces',f,'Vertices',v,'FaceColor',pc_thal_colors(i,:),'Clipping','off','EdgeColor','none');
end
for i = 1:921
direction = [0 0 1];
rotate(pc_patches(i),direction,46)

% direction = [1 0 0];
% rotate(pc_patches(i),direction,13)
end

for i = 1:921

direction = [1 0 0];
rotate(pc_patches(i),direction,-13)
end