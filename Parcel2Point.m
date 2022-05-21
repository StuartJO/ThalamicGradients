function Parcel2Point(surface,roi_color,Surfdata,EndLocation)

[BOUNDARY,BOUNDARY_ROI_ID] = findROIboundaries(surface.vertices,lh_faces,lh_rand500,'midpoint');

Nrois = length(Surfdata);

C1 = rescale(coeff{3}(1:250,1),ylimits(1),ylimits(2))*-1;
C2 = rescale(coeff{3}(1:250,2),zlimits(1),zlimits(2));

for i =1:Nrois
   [x{i},y{i}] = circle(C1(i),C2(i),2,length(BOUNDARY{i+1})); 
end


for i = 1:Nrois
bnew= find_point_on_line(BOUNDARY{i+1},[zeros(length(x{i}),1) x{i} y{i}],0);
hold on
p3(i)=fill3(bnew(:,1),bnew(:,2),bnew(:,3),roi_color(i,:),'FaceLighting','gouraud','LineWidth',2,'Clipping','off');
material dull

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
for i = 1:Nrois
  bnew= find_point_on_line(BOUNDARY{i+1},[zeros(length(x{i}),1) x{i} y{i}],r(j));  
    
p3(i).Vertices = bnew;
ylim(ylimits)
zlim(zlimits)
end
%print(['TestParcel2Point',num2str(j),'.png'],'-dpng')
pause(.1)
end
