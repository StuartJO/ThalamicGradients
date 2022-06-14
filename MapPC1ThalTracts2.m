figure('Position',[0 0 1920 1080])
addpath F:\Documents\GitHub\plotSurfaceROIBoundary
axes('Position',[0.0177 0.0241 0.9698 0.9657])

load('fsaverage_surface_data.mat')
%plotSurfaceROIBoundary(surface,lh_rand500,IT,'midpoint',IT_colormap,2);
lh_inflated_verts_med = lh_inflated_verts;

theta=3;
lh_inflated_verts_med(:,1) = lh_inflated_verts(:,1)*cos(theta) - lh_inflated_verts(:,2)*sin(theta);
lh_inflated_verts_med(:,2) = lh_inflated_verts(:,1)*sin(theta) + lh_inflated_verts(:,2)*cos(theta) - 220;
lh_inflated_verts_med(:,3) = lh_inflated_verts(:,3);

BOUNDARY = findROIboundaries(lh_inflated_verts,lh_faces,lh_rand500,'midpoint');



for i = 1:251
hold on
pIT(i)=fill3(BOUNDARY{i}(:,1),BOUNDARY{i}(:,2),BOUNDARY{i}(:,3),[.5 .5 .5],'FaceLighting','gouraud','LineWidth',2,'Clipping','off');
material dull

end

BOUNDARY2 = findROIboundaries(lh_inflated_verts_med,lh_faces,lh_rand500,'midpoint');
for i = 1:251
hold on
pIT2(i)=fill3(BOUNDARY2{i}(:,1),BOUNDARY2{i}(:,2),BOUNDARY2{i}(:,3),[.5 .5 .5],'FaceLighting','gouraud','LineWidth',2,'Clipping','off');
material dull

end


view([-90 0])
axis off
axis image

ylimits = ylim;
zlimits = zlim;
xlimits = xlim;


[~,score_sort] = sort(zscore(score{3}(:,1)),'descend');

cort_data = norm(score_sort(1),:);

tract_cmap = parula(256);

 tract_color_ind = ceil(rescale(cort_data,1,size(tract_cmap,1)));
 tract_color = tract_cmap(tract_color_ind,:);
 tract_color = [.5 .5 .5; tract_color];
 
Points2Sample = round(linspace(1,921,120));

for j = 1:120

cort_data = norm(score_sort(Points2Sample(j)),:);

tract_cmap = parula(256);

 tract_color_ind = ceil(rescale(cort_data,1,size(tract_cmap,1)));
 tract_color = tract_cmap(tract_color_ind,:);
 tract_color = [.5 .5 .5; tract_color];
for i = 1:251
pIT(i).FaceColor = tract_color(i,:);
pIT2(i).FaceColor = tract_color(i,:);
end

print(['./NewGIFs/PC1_tracts_',num2str(j),'.png'],'-dpng')

end