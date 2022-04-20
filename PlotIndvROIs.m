



surface.vertices = lh_inflated_verts;
surface.faces = lh_faces;
plotSurfaceROIBoundary(surface,lh_rand500,1:250,'midpoint',turbo(250),1,2);

% The following options set up the patch object to look pretty. This works
% well for the left hemisphere (medial and lateral). Change the inputs to 
% 'view' to see the brain from different angles ([-90 0] for left and [90 0]
% for right I find works well)

camlight(80,-10);
camlight(-80,-10);
view([-90 0])
axis off
axis tight
axis equal

for i = 1:250
    roi_verts_ind = find(lh_rand500==i);
    roi_verts = surface.vertices(roi_verts_ind,:);
    roi_faces_ind = sum(ismember(lh_faces,roi_verts_ind),2)==3;
    roi_faces = lh_faces(roi_faces_ind,:);
    roi_faces_new = changem(roi_faces,1:length(roi_verts_ind),roi_verts_ind);
    verts_rois{i} = roi_verts;
    faces_rois{i} = roi_faces_new;
end

[BOUNDARY,BOUNDARY_ROI_ID] = findROIboundaries(surface.vertices,lh_faces,lh_rand500,'midpoint');

cmap = turbo(250);
figure
for i = 1:250
    
   p(i) = patch('Vertices',verts_rois{i},'Faces',faces_rois{i});
set(p(i),'EdgeColor','none','FaceColor',cmap(i,:),'Clipping','off');
p(i).FaceLighting = 'gouraud';
material dull
hold on

boundary_plot.boundary(i) = plot3(BOUNDARY{i}(:,1), BOUNDARY{i}(:,2), BOUNDARY{i}(:,3), 'Color', 'k', 'LineWidth',2,'Clipping','off');
end
camlight(80,-10);
camlight(-80,-10);
view([-90 0])
axis off
axis tight
axis equal

ylimits = ylim;
zlimits = zlim;

C1 = rescale(coeff{3}(1:250,1),ylimits(1),ylimits(2))*-1;
C2 = rescale(coeff{3}(1:250,2),zlimits(1),zlimits(2));

for r = 0:.01:1
%r = .8;
for i = 1:250
    
   v = verts_rois{i};
   newPos = repmat([0,C1(i),C2(i)],length(v),1);
   vnew = find_point_on_line(v,newPos,r);
    newPos2 = repmat([0,C1(i),C2(i)],length(BOUNDARY{i+1}),1);
    bnew= find_point_on_line(BOUNDARY{i+1},newPos2,r);
    
p(i).Vertices = vnew;
boundary_plot.boundary(i).XData = bnew(:,1);
boundary_plot.boundary(i).YData = bnew(:,2);
boundary_plot.boundary(i).ZData = bnew(:,3);
ylim(ylimits)
zlim(zlimits)
end
pause(.1)
end

%%%

for i =1:250
   [x{i},y{i}] = circle(C1(i),C2(i),2,length(BOUNDARY{i+1})); 
end

figure
for i = 1:250
%newPos2 = repmat([0,C1(i),C2(i)],length(BOUNDARY{i+1}),1);
%bnew= find_point_on_line(BOUNDARY{i+1},newPos2,r);
bnew= find_point_on_line(BOUNDARY{i+1},[zeros(length(x{i}),1) x{i} y{i}],0);
hold on
p3(i)=fill3(bnew(:,1),bnew(:,2),bnew(:,3),cmap(i,:),'FaceLighting','gouraud','LineWidth',2,'Clipping','off');
material dull
    
% B.boundary(i)= plot3(bnew(:,1), bnew(:,2), bnew(:,3), 'Color', 'k', 'LineWidth',2,'Clipping','off');

% boundary_plot.boundary(i).XData = bnew(:,1);
% boundary_plot.boundary(i).YData = bnew(:,2);
% boundary_plot.boundary(i).ZData = bnew(:,3);
end
camlight(80,-10);
camlight(-80,-10);
view([-90 0])
axis off
axis tight
axis equal

ylimits = ylim;
zlimits = zlim;

r =0:.01:1;

for j = 1:length(r)
%r = .8;
for i = 1:250
  bnew= find_point_on_line(BOUNDARY{i+1},[zeros(length(x{i}),1) x{i} y{i}],r(j));  
    
p3(i).Vertices = bnew;
ylim(ylimits)
zlim(zlimits)
end
print(['TestParcel2Point',num2str(j),'.png'],'-dpng')
pause(.1)
end

% 
% 
% figure
% scatter(coeff{3}(1:250,1),coeff{3}(1:250,2),10,1:250,'filled')
% 
% 
% p2 = patch('Vertices',lh_verts,'Faces',lh_faces);
% set(p2,'EdgeColor','none','FaceColor',[.5 .5 .5],'Clipping','off','FaceAlpha',.2,'FaceLighting','gouraud');
