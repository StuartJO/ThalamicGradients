figure('Position',[0 0 1920 1080])

% IT = neuromap_parc(:,24);
% IT_colormap = cbrewer('div','RdBu',256);

[~,score_sort] = sort(score{3}(:,1),'descend');

cort_data = norm(score_sort(1),:);

tract_cmap = parula(256);

 tract_color_ind = ceil(rescale(cort_data,1,size(tract_cmap,1)));
 tract_color = tract_cmap(tract_color_ind,:);
 tract_color = [.5 .5 .5; tract_color];
 

surface.vertices = lh_inflated_verts;
surface.faces = lh_faces;
%plotSurfaceROIBoundary(surface,lh_rand500,IT,'midpoint',IT_colormap,2);

BOUNDARY = findROIboundaries(lh_inflated_verts,lh_faces,lh_rand500,'midpoint');
for i = 1:251
hold on
pIT(i)=fill3(BOUNDARY{i}(:,1),BOUNDARY{i}(:,2),BOUNDARY{i}(:,3),tract_color(i,:),'FaceLighting','gouraud','LineWidth',2,'Clipping','off');
material dull

end

camlight(80,-10);
camlight(-80,-10);
view([-90 0])
axis off
axis tight
axis equal

Points2Sample = round(linspace(1,921,120));

for j = 1:120

cort_data = norm(score_sort(Points2Sample(j)),:);

tract_cmap = parula(256);

 tract_color_ind = ceil(rescale(cort_data,1,size(tract_cmap,1)));
 tract_color = tract_cmap(tract_color_ind,:);
 tract_color = [.5 .5 .5; tract_color];
for i = 1:251
pIT(i).FaceColor = tract_color(i,:);
end

%print(['./GIF/MAPPC1_tracts_',num2str(j),'.png'],'-dpng')
pause(.1)
end

view([90 0])

for j = 1:120

cort_data = norm(Points2Sample(j),:);

tract_cmap = parula(256);

 tract_color_ind = ceil(rescale(cort_data,1,size(tract_cmap,1)));
 tract_color = tract_cmap(tract_color_ind,:);
 tract_color = [.5 .5 .5; tract_color];
for i = 1:251
pIT(i).FaceColor = tract_color(i,:);
end

%print(['./GIF/MAPPC1_tracts_medial_',num2str(j),'.png'],'-dpng')
pause(.1)
end