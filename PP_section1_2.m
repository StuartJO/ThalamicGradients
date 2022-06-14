load('Sub76_ThalData.mat')

load('fsaverage_surface_data.mat')
load('lh.avgMapping_allSub_RF_ANTs_MNI152_orig_to_fsaverage.mat');
lh_mni_vox = MNI_mm2vox(ras','mm');
load('rh.avgMapping_allSub_RF_ANTs_MNI152_orig_to_fsaverage.mat');
rh_mni_vox = MNI_mm2vox(ras','mm');
MNIsurface.vertices = lh_mni_vox;
MNIsurface.faces = lh_faces;


figure('Position',[0 0 1920 1080])
MNIsurface.vertices = lh_mni_vox;
MNIsurface.faces = lh_faces;
%pMNIsurface = plotSurfaceROIBoundary(MNIsurface,lh_rand500,1:250,'none',turbo(250),2);

pMNIsurface = patch(MNIsurface);
set(pMNIsurface,'EdgeColor','none','FaceColor',[.5 .5 .5],'Clipping','off');
pMNIsurface.FaceLighting = 'gouraud';
material dull

camlight(80,-10);
camlight(-80,-10);
view([90 0])

axis off
axis tight
axis equal
axis vis3d

hold on

%TR = stlread_('LeftThalMeshSmoothManual.stl');

TR = stlread('LeftThalMeshSmoothManual.stl');

Thalsurf.faces = TR.ConnectivityList;
Thalsurf.vertices = MNI_mm2vox(TR.Points,'mm');

% [v1, f1] = subdivide_tri(TR.Points, TR.ConnectivityList);

% Thalsurf.faces = f1;
% Thalsurf.vertices = MNI_mm2vox(v1,'mm');

pThalsurf = patch(Thalsurf);
set(pThalsurf,'EdgeColor','none','FaceColor',[1 .647 0],'Clipping','off');
pThalsurf.FaceLighting = 'gouraud';
material dull

xlimits = xlim;
ylimits = ylim;
zlimits = zlim;
alphavals = linspace(1,.025,30);
i = 30;
pMNIsurface.FaceAlpha = alphavals(i);

end_xlimits = [133.3088  143.8052];
end_ylimits = [133.1712  162.0478];
end_zlimits = [81.3852  105.1673];


t= 60;
views = [linspace(90,136,t)' linspace(0,13,t)'];
thal_alphavals = linspace(1,.025,t);
xlimits_range = [linspace(xlimits(1),end_xlimits(1),t)' linspace(xlimits(2),end_xlimits(2),t)'];
ylimits_range = [linspace(ylimits(1),end_ylimits(1),t)' linspace(ylimits(2),end_ylimits(2),t)'];
zlimits_range = [linspace(zlimits(1),end_zlimits(1),t)' linspace(zlimits(2),end_zlimits(2),t)'];
for i = 1:t
    pThalsurf.FaceAlpha = 1;
    view(views(i,:))
    xlim(xlimits_range(i,:));
    ylim(ylimits_range(i,:));
    zlim(zlimits_range(i,:));
    pause(.1)
    print(['./GIF/PP1_S3_',num2str(i),'.png'],'-dpng')
end


[~,ThalTianParc] = read_nifti('TianS4_LeftThal.nii');
[I1,I2,I3] = ind2sub(size(ThalTianParc),find(ThalTianParc~=0));
ThalVerts = round(Thalsurf.vertices);
ThalTianParcParcID=ThalTianParc(find(ThalTianParc~=0));
[~,ThalVerts_parcvoxIND] = min(pdist2(ThalVerts,[I1,I2,I3]),[],2);
Tian = ThalTianParcParcID(ThalVerts_parcvoxIND);

[~,ThalFSLParc] = read_nifti('striatum6ANDthalamus14_config.nii');
[I1,I2,I3] = ind2sub(size(ThalFSLParc),find(ThalFSLParc>3));
ThalVerts = round(Thalsurf.vertices);
ThalFSLParcParcID=ThalFSLParc(find(ThalFSLParc>3));
[~,ThalVerts_parcvoxIND] = min(pdist2(ThalVerts,[I1,I2,I3]),[],2);
FSL = ThalFSLParcParcID(ThalVerts_parcvoxIND);
FSL = changem(FSL,1:length(unique(FSL)),unique(FSL));

[~,ThalMorelParc] = read_nifti('morel_lh.nii');
ThalMorelParcParcID=ThalMorelParc(find(ThalMorelParc~=0));

ThalMorelParcParcID_unique = unique(ThalMorelParcParcID);
counts = histc(ThalMorelParcParcID(:), ThalMorelParcParcID_unique);

Parcs2Use = ThalMorelParcParcID_unique(counts>150);

ThalMorelParcParcID2=ThalMorelParc(find(ismember(ThalMorelParc,Parcs2Use)));
[I1,I2,I3] = ind2sub(size(ThalMorelParc),find(ismember(ThalMorelParc,Parcs2Use)));
ThalVerts = round(Thalsurf.vertices);


[~,ThalVerts_parcvoxIND] = min(pdist2(ThalVerts,[I1,I2,I3]),[],2);
ThalVertID = ThalMorelParcParcID2(ThalVerts_parcvoxIND);
MOREL = changem(ThalVertID,1:length(unique(ThalVertID)),unique(ThalVertID));

%PARC = [ones(length(Thalsurf.vertices),1) Tian FSL MOREL];

PARC = [ones(length(Thalsurf.vertices),1) MOREL FSL Tian MOREL];

for i = 1:size(PARC,2)
    Nrois = max(PARC(:,i));
    parc_dist_temp = zeros(length(PARC),1);
    for j = 1:Nrois
        parc_inds = PARC(:,i)==j;
        parc_verts = Thalsurf.vertices(parc_inds,:);
        parc_Dists = squareform(pdist(parc_verts));
        mean_parc_Dists = mean(parc_Dists);
        [~,closest_vert] = min(mean_parc_Dists);
        dist2centre = parc_Dists(:,closest_vert);
       parc_dist_temp(parc_inds) = dist2centre;
    end
    parc_dist{i} = parc_dist_temp;
    parc_dist_norm{i} = parc_dist_temp./max(parc_dist_temp);
end
hot_cmap = hot(9);
cmap = [[1 .647 0];lines(15); parula(9);hot_cmap(2:9,:);lines(15);];

cmap_length = size(cmap,1);

PARCnew = PARC(:,1);
for i = 2:size(PARC,2)
PARCnew(:,i) = PARC(:,i)+max(PARCnew(:,i-1));
end

parc_dist2centre = zeros(size(PARCnew,1),4);
for i = 1:5
parc_dist2centre(:,i) = parc_dist{i};
end

PARCnew(PARCnew(:,1)==0,:) = 0;

parc = PARCnew(:,1);

pThalsurf.FaceColor = 'flat';

COUNT = 0;

MOREL_BOUN_FACES = findROIboundaries(Thalsurf.vertices,Thalsurf.faces,MOREL,'faces');

for i = [2:5]
   
    maxdist = max(parc_dist2centre(:,i));
    
    distThr = linspace(0,maxdist,60);
    
    for j = 1:60
        addnewparc_ind = parc_dist2centre(:,i)<distThr(j);
        parc(addnewparc_ind) = PARCnew(addnewparc_ind,i);
        if i == 2
       FaceVertexCData = makeFaceVertexCData(Thalsurf.vertices,Thalsurf.faces,parc,parc,cmap,[1 cmap_length],1);
        else
             
      FaceVertexCData = makeFaceVertexCData(Thalsurf.vertices,Thalsurf.faces,parc,parc,cmap,[1 cmap_length],1,[.5 .5 .5],[.75 .75 .75]);
      FaceVertexCData(MOREL_BOUN_FACES,1)=0;
      FaceVertexCData(MOREL_BOUN_FACES,2)=0;
      FaceVertexCData(MOREL_BOUN_FACES,3)=0;
        end
    pThalsurf.FaceVertexCData = FaceVertexCData;

        pause(.1)
        print(['./GIF/PP1_S',num2str(i+2),'_',num2str(j),'.png'],'-dpng')
    end
    
    
end
% 
% 
% fc1 = gifti('hcp.embed.grad_1.L.fsa.func.gii');
% fc1_data = fc1.cdata;
% 
% FaceVertexCData = makeFaceVertexCData(MNIsurface.vertices,MNIsurface.faces,lh_rand500,fc1_data,viridisplus(256));
% 
% r = linspace(0,1,60);
% 
% FaceVertexCData_gray = FaceVertexCData;
% FaceVertexCData_gray(:) = .5;
% alphavals = linspace(.025,1,60);
% 
% for i = 1:60
% NewFaceColors = find_point_on_line(FaceVertexCData_gray,FaceVertexCData,r(i));
%     view(views(61-i,:))
%     xlim(xlimits_range(61-i,:));
%     ylim(ylimits_range(61-i,:));
%     zlim(zlimits_range(61-i,:));
%     pMNIsurface.FaceAlpha = alphavals(i);
%     pMNIsurface.FaceVertexCData = NewFaceColors;
%     pMNIsurface.FaceColor='interp';
%     pause(.1)
%     print(['./GIF/PP1_S7_',num2str(i),'.png'],'-dpng')
% end
