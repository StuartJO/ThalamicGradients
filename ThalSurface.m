%ThalSurf = gifti('lh.LeftThalMeshSmooth2.surf.gii');

TR = stlread('LeftThalMeshSmoothManual.stl');

%Thalsurf.vertices = ThalSurf.vertices;
%Thalsurf.faces = ThalSurf.faces;
Thalsurf.faces = TR.ConnectivityList;
Thalsurf.vertices = MNI_mm2vox(TR.Points,'mm');
%Thalsurf.vertices = MNI_mm2vox(ThalSurf.vertices,'mm');

surface2.vertices = lh_mni_vox;
[p2,b,BOUNDARY] = plotSurfaceROIBoundary(surface2,lh_rand500,1:250,'midpoint',turbo(250),2);

camlight(80,-10);
camlight(-80,-10);
view([90 0])
xlimits = xlim;
ylimits = ylim;
zlimits = zlim;
%axis off
axis tight
axis equal
axis vis3d
hold on

p2.FaceAlpha=.1;


thalpatch = patch(Thalsurf);
set(thalpatch,'EdgeColor','none','FaceColor',[.5 .5 .5],'Clipping','off');
thalpatch.FaceLighting = 'gouraud';
material dull
camlight(80,-10);
camlight(-80,-10);
view([90 0])
axis tight
axis equal
axis vis3d

[~,ThalTianParc] = read_nifti('TianS4_LeftThal.nii');
[I1,I2,I3] = ind2sub(size(ThalTianParc),find(ThalTianParc~=0));
ThalVerts = round(Thalsurf.vertices);
ThalTianParcParcID=ThalTianParc(find(ThalTianParc~=0));
[~,ThalVerts_parcvoxIND] = min(pdist2(ThalVerts,[I1,I2,I3]),[],2);
ThalVertID = ThalTianParcParcID(ThalVerts_parcvoxIND);

[thalpatch] = plotSurfaceROIBoundary(Thalsurf,ThalVertID,1:8,'midpoint',turbo(8),2);
camlight(80,-10);
camlight(-80,-10);
view([136 13])
axis tight
axis equal
axis vis3d



[~,ThalMorelParc] = read_nifti('morel_lh.nii');
ThalMorelParcParcID=ThalMorelParc(find(ThalMorelParc~=0));

ThalMorelParcParcID_unique = unique(ThalMorelParcParcID);
counts = histc(ThalMorelParcParcID(:), ThalMorelParcParcID_unique);

Parcs2Use = ThalMorelParcParcID_unique(counts>100);

ThalMorelParcParcID2=ThalMorelParc(find(ismember(ThalMorelParc,Parcs2Use)));
[I1,I2,I3] = ind2sub(size(ThalMorelParc),find(ismember(ThalMorelParc,Parcs2Use)));
ThalVerts = round(Thalsurf.vertices);


[~,ThalVerts_parcvoxIND] = min(pdist2(ThalVerts,[I1,I2,I3]),[],2);
ThalVertID = ThalMorelParcParcID2(ThalVerts_parcvoxIND);
ThalVertIDNew = changem(ThalVertID,unique(ThalVertID),1:length(unique(ThalVertID)));

[thalpatch] = plotSurfaceROIBoundary(Thalsurf,ThalVertIDNew,ThalVertIDNew,'midpoint',turbo(20),2);
camlight(80,-10);
camlight(-80,-10);
view([136 13])
axis tight
axis equal
axis vis3d


[~,ThalFSLParc] = read_nifti('striatum6ANDthalamus14_config.nii');
[I1,I2,I3] = ind2sub(size(ThalFSLParc),find(ThalFSLParc>3));
ThalVerts = round(Thalsurf.vertices);
ThalFSLParcParcID=ThalFSLParc(find(ThalFSLParc>3));
[~,ThalVerts_parcvoxIND] = min(pdist2(ThalVerts,[I1,I2,I3]),[],2);
ThalVertID = ThalFSLParcParcID(ThalVerts_parcvoxIND);

[thalpatch] = plotSurfaceROIBoundary(Thalsurf,ThalVertID,ThalVertID,'midpoint',lines(20),2);
camlight(80,-10);
camlight(-80,-10);
view([136 13])
axis tight
axis equal
axis vis3d