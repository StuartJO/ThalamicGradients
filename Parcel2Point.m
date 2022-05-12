function Parcel2Point(surface,vertex_id,Surfdata,EndLocation)

[BOUNDARY,BOUNDARY_ROI_ID] = findROIboundaries(surface.vertices,lh_faces,lh_rand500,'midpoint');

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
