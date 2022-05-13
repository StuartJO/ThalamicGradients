load('fsaverage_surface_data.mat')
load('lh.avgMapping_allSub_RF_ANTs_MNI152_orig_to_fsaverage.mat');
lh_mni_vox = MNI_mm2vox(ras','mm');
load('rh.avgMapping_allSub_RF_ANTs_MNI152_orig_to_fsaverage.mat');
rh_mni_vox = MNI_mm2vox(ras','mm');
MNIsurface.vertices = lh_mni_vox;
MNIsurface.faces = lh_faces;

MNIsurfaceR.vertices = rh_mni_vox;
MNIsurfaceR.faces = rh_faces;

tracks = read_mrtrix_tracks ('192035inMNI_thal_seed_620.tck');

pMNIsurface = patch(MNIsurface);
set(pMNIsurface,'EdgeColor','none','FaceColor',[.5 .5 .5],'Clipping','off');
pMNIsurface.FaceLighting = 'gouraud';
material dull

view([180 0])

axis off
axis tight
axis equal
axis vis3d

hold on
pMNIsurfaceR = patch(MNIsurfaceR);
set(pMNIsurfaceR,'EdgeColor','none','FaceColor',[.5 .5 .5],'Clipping','off');
pMNIsurfaceR.FaceLighting = 'gouraud';
material dull

clight = camlight('headlight');

xlimits = xlim;
ylimits = ylim;
zlimits = zlim;

end_xlimits = [133.3088  143.8052];
end_ylimits = [133.1712  162.0478];
end_zlimits = [81.3852  105.1673];

t= 60;
views = [linspace(180,136,t)' linspace(0,13,t)'];
xlimits_range = [linspace(xlimits(1),end_xlimits(1),t)' linspace(xlimits(2),end_xlimits(2),t)'];
ylimits_range = [linspace(ylimits(1),end_ylimits(1),t)' linspace(ylimits(2),end_ylimits(2),t)'];
zlimits_range = [linspace(zlimits(1),end_zlimits(1),t)' linspace(zlimits(2),end_zlimits(2),t)'];

TR = stlread('LeftThalMeshSmoothManual.stl');
Thalsurf.faces = TR.ConnectivityList;
Thalsurf.vertices = MNI_mm2vox(TR.Points,'mm');

pThalsurf = patch(Thalsurf);
set(pThalsurf,'EdgeColor','none','FaceColor',[1 .647 0],'Clipping','off');
pThalsurf.FaceLighting = 'gouraud';
material dull


Cortalphavals = [1 1 1 1 1 1 1 1 1 1 linspace(1,.025,50)];
CortalphavalsR = [1 1 1 1 1 1 1 1 1 1 linspace(1,0,50)];
Thalalphavals = [0 0 0 0 0 0 0 0 0 0 ones(1,50)];

for i = 1:60
    pThalsurf.FaceAlpha = Thalalphavals(i);
    pMNIsurface.FaceAlpha = Cortalphavals(i);
    pMNIsurfaceR.FaceAlpha = Cortalphavals(i);
    view(views(i,:))
    xlim(xlimits_range(i,:));
    ylim(ylimits_range(i,:));
    zlim(zlimits_range(i,:));
    pause(.1)    
    delete(clight)
    clight = camlight('headlight');
    pause(.1)
end
% camlight(80,-10);
% camlight(-80,-10);
%%

pThalsurf.FaceAlpha=.025;
[x,y,z]=sphere(64);
sphere_mesh = surf2patch(x,y,z,'triangles');
load('MNI_Seed_voxelData.mat')
seed_voxel_coords = seeds_vox(logical(seed_ind),:);
for i = 1:921
    sphere_mesh_roi = sphere_mesh;
    sphere_mesh_roi.vertices = (sphere_mesh.vertices.*.15) + seed_voxel_coords(i,:);
    sSeeds(i) = patch(sphere_mesh_roi,'EdgeColor','none','FaceColor',[0 0 0]./255,'Clipping','off','FaceLighting','gouraud');
    material dull
end





z_loc = 82;
brain_slice = brain(:,:,z_loc)';
colormap(gray(256))
slicesurf = surf(repmat(z_loc, size(brain_slice)), brain_slice);
slicesurf.EdgeAlpha=0;
slicesurf.FaceAlpha=0;
slicesurf.FaceLighting='none';
slicesurf.Clipping='off';

xlim(end_xlimits)
ylim(end_ylimits)
zlim(end_zlimits)
%axis tight


seed_nearest = 83.25;
seeds2plot_ind = (seed_voxel_coords(:,3)==seed_nearest);


slice_xlimits = [133.3088  143.8052];
slice_ylimits = [92.5534  126.3029];
slice_zlimits = [116.0609  143.8561];

t= 60;
views = [linspace(136,0,t)' linspace(13,90,t)'];
xlimits_range = [linspace(end_xlimits(1),slice_xlimits(1),t)' linspace(end_xlimits(2),slice_xlimits(2),t)'];
ylimits_range = [linspace(end_ylimits(1),slice_ylimits(1),t)' linspace(end_ylimits(2),slice_ylimits(2),t)'];
zlimits_range = [linspace(end_zlimits(1),slice_zlimits(1),t)' linspace(end_zlimits(2),slice_zlimits(2),t)'];

    pThalsurf.FaceAlpha = 0;
    pMNIsurface.FaceAlpha = 0;
    pMNIsurfaceR.FaceAlpha = 0;

    
voxel_seeds = seeds_vox(logical(seed_ind),:);
thr=2.1;
GRAD_DATA = pc_thal;
inmask = find(thalmask==1);
[X,Y,Z] = ind2sub(size(thalmask),inmask);
ThalMaskDists = pdist2([X Y Z],voxel_seeds);
[ThalMask_nearest_seedmm,ThalMask_nearest_seed] = min(ThalMaskDists');
ThalVoxClust = changem(ThalMask_nearest_seed,GRAD_DATA,1:length(voxel_seeds));
thalmaskclust = thalmask;
thalmaskclust(thalmaskclust==0) = NaN;
thalmaskclust(inmask) = ThalVoxClust;
thalmaskclust(inmask(ThalMask_nearest_seedmm>thr)) = NaN;

thalmaskclust_scaled = rescale(thalmaskclust,257,512);

thalslice = thalmaskclust_scaled(:,:,z_loc)';

newCmap = [gray(256); turbo(256)];

newslicevals1 = rescale(brain_slice,1,256);

newslicevals1(find(~isnan(newslicevals2))) = thalslice(find(~isnan(newslicevals2)));

colormap(newCmap)
slicesurf.CData=newslicevals1;

gene_data = data_type{2};
zscore(gene_data(:,2))






for i = 1:5000
tracts_MNI{i} = MNI_mm2vox(tracks.data{i},'mm');

%inside = PointInsideVolume(tracts_MNI{i}, lh_faces, lh_mni_vox);

end



for i = 1:5000%length(OK_STREAMLINES)
    S = tracts_MNI{i};
    
    if size(S,1) > 1
    
    [~,SDIR] = max(sum(abs(diff(S))));    
    
    if SDIR == 1
        SCOL = [1 0 0 .2];
        col = 'r';
    elseif SDIR == 2
        SCOL = [0 1 0 .2];
        col = 'g';
    else
        SCOL = 	[0 0 1 .2];
        col = 'b';
    end
    
    
    x = S(:,1);
    y = S(:,2); 
    z = S(:,3); 

    %handles.handle_plotCD(i) = plot3(x,y,z,'LineWidth',1,'Color',SCOL,'Clipping','off');
    %handles.handle_plotCD(i) = plot3t(x,y,z,.5,col);
    StepRate = 1;

    handles.handle_plotCD(i) = patch([S(:,1)' NaN],[S(:,2)' NaN],[S(:,3)' NaN],0); 
    joint = ([S(2:end,:); NaN NaN NaN]-S);
    joint(end,:) = joint(end-1,:);
    temp_joint = joint;
    joint(:,1) = temp_joint(:,1);
    joint(:,2) = temp_joint(:,2);
    cdata = [abs(joint./StepRate); NaN NaN NaN];
    cdata = reshape(cdata,length(cdata),1,3);
    set(handles.handle_plotCD(i),'CData', cdata, 'EdgeColor','interp','FaceColor','interp','Clipping','off') 
    
    hold on
    
    end
end