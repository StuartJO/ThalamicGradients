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
for i = 1:30
pMNIsurface.FaceAlpha = alphavals(i);
print(['./GIF/PP1_S1_',num2str(i),'.png'],'-dpng')
end



%tracks = read_mrtrix_tracks ('combined_tracts_reduced_MNI.tck');
tracks = read_mrtrix_tracks ('192035inMNI_RAND_THAL_TRACTS2.tck');

for i = 1:length(tracks.data)
tracts_MNI{i} = MNI_mm2vox(tracks.data{i},'mm');

%inside = PointInsideVolume(tracts_MNI{i}, lh_faces, lh_mni_vox);

end
t = 60;
for j = 1:t
    if exist('streamline_handle','var')
    delete(streamline_handle)
    end
    for i = 1:length(tracts_MNI)
        s = tracts_MNI{i};
        streamline_length = size(s,1);

        Sr = round(find_point_on_line(1,streamline_length,j/60));
        
        S = s(1:Sr,:);
        
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


%         x = S(:,1);
%         y = S(:,2); 
%         z = S(:,3); 

        %handles.handle_plotCD(i) = plot3(x,y,z,'LineWidth',1,'Color',SCOL,'Clipping','off');
        %handles.handle_plotCD(i) = plot3t(x,y,z,.5,col);
        StepRate = 1;

        streamline_handle(i) = patch([S(:,1)' NaN],[S(:,2)' NaN],[S(:,3)' NaN],0); 
        joint = ([S(2:end,:); NaN NaN NaN]-S);
        joint(end,:) = joint(end-1,:);
        temp_joint = joint;
        joint(:,1) = temp_joint(:,1);
        joint(:,2) = temp_joint(:,2);
        cdata = [abs(joint./StepRate); NaN NaN NaN];
        cdata = reshape(cdata,length(cdata),1,3);
        set(streamline_handle(i),'CData', cdata, 'EdgeColor','interp','FaceColor','interp','Clipping','off') 

        hold on

        end
    end
xlim(xlimits);
ylim(ylimits);
zlim(zlimits);
    pause(.1)
    print(['./GIF/PP1_S2_',num2str(j),'.png'],'-dpng')
end



print(['./GIF/PP1_S2.png'],'-dpng')

delete(streamline_handle)

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


delete(pThalsurf)

[~,ThalTianParc] = read_nifti('TianS4_LeftThal.nii');
[I1,I2,I3] = ind2sub(size(ThalTianParc),find(ThalTianParc~=0));
ThalVerts = round(Thalsurf.vertices);
ThalTianParcParcID=ThalTianParc(find(ThalTianParc~=0));
[~,ThalVerts_parcvoxIND] = min(pdist2(ThalVerts,[I1,I2,I3]),[],2);
ThalVertID = ThalTianParcParcID(ThalVerts_parcvoxIND);

[thalpatch,thalpatchboundary] = plotSurfaceROIBoundary(Thalsurf,ThalVertID,ThalVertID,'faces',turbo(8),10);

print(['./GIF/PP1_S4.png'],'-dpng')

delete(thalpatch)
%delete(thalpatchboundary.boundary)

[~,ThalFSLParc] = read_nifti('striatum6ANDthalamus14_config.nii');
[I1,I2,I3] = ind2sub(size(ThalFSLParc),find(ThalFSLParc>3));
ThalVerts = round(Thalsurf.vertices);
ThalFSLParcParcID=ThalFSLParc(find(ThalFSLParc>3));
[~,ThalVerts_parcvoxIND] = min(pdist2(ThalVerts,[I1,I2,I3]),[],2);
ThalVertID = ThalFSLParcParcID(ThalVerts_parcvoxIND);

[thalpatch,thalpatchboundary] = plotSurfaceROIBoundary(Thalsurf,ThalVertID,ThalVertID,'faces',lines(20),10);

print(['./GIF/PP1_S5.png'],'-dpng')

delete(thalpatch)
%delete(thalpatchboundary.boundary)

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
ThalVertIDNew = changem(ThalVertID,unique(ThalVertID),1:length(unique(ThalVertID)));

[thalpatch,thalpatchboundary] = plotSurfaceROIBoundary(Thalsurf,ThalVertIDNew,ThalVertIDNew,'midpoint',lines(20),10);

print(['./GIF/PP1_S6.png'],'-dpng')

fc1 = gifti('hcp.embed.grad_1.L.fsa.func.gii');
fc1_data = fc1.cdata;

FaceVertexCData = makeFaceVertexCData(MNIsurface.vertices,MNIsurface.faces,lh_rand500,fc1_data,viridisplus(256));

r = linspace(0,1,60);

FaceVertexCData_gray = FaceVertexCData;
FaceVertexCData_gray(:) = .5;
alphavals = linspace(.025,1,60);

for i = 1:60
NewFaceColors = find_point_on_line(FaceVertexCData_gray,FaceVertexCData,r(i));
    view(views(61-i,:))
    xlim(xlimits_range(61-i,:));
    ylim(ylimits_range(61-i,:));
    zlim(zlimits_range(61-i,:));
    pMNIsurface.FaceAlpha = alphavals(i);
    pMNIsurface.FaceVertexCData = NewFaceColors;
    pMNIsurface.FaceColor='interp';
    pause(.1)
    print(['./GIF/PP1_S7_',num2str(i),'.png'],'-dpng')
end
