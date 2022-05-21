figure('Position',[0 0 1920 1080])

IT = neuromap_parc(:,24);
IT_colormap = cbrewer('div','RdBu',256);



 IT_color_ind = ceil(rescale(IT,1,size(IT_colormap,1)));
 IT_color = IT_colormap(IT_color_ind,:);
 IT_color = [.5 .5 .5; IT_color];
 

surface.vertices = lh_inflated_verts;
surface.faces = lh_faces;
%plotSurfaceROIBoundary(surface,lh_rand500,IT,'midpoint',IT_colormap,2);

BOUNDARY = findROIboundaries(lh_inflated_verts,lh_faces,lh_rand500,'midpoint');
for i = 1:251
hold on
pIT(i)=fill3(BOUNDARY{i}(:,1),BOUNDARY{i}(:,2),BOUNDARY{i}(:,3),IT_color(i,:),'FaceLighting','gouraud','LineWidth',2,'Clipping','off');
material dull

end

for i = 1:251
mean_expanded = mean(BOUNDARY{i}*1.1)-mean(BOUNDARY{i});

BOUNDARY_exploded{i} = BOUNDARY{i}+mean_expanded;
end

camlight(80,-10);
camlight(-80,-10);
view([-90 0])
axis off
axis tight
axis equal

ylimits = ylim;
zlimits = zlim;
xlimits = xlim;

for i = 1:251
hold on
pIT(i).FaceColor=[.5 .5 .5];
end

print(['./GIF/PPEB_S1.png'],'-dpng')

for i = 1:251
hold on
pIT(i).FaceColor=IT_color(i,:);
end

print(['./GIF/PPEB_S2.png'],'-dpng')
 PC1_color = turbo(256);
 PC1_cort_color = value2color(coeff{3}(1:250,1),PC1_color);
 
 PC1_cort_color = [.5 .5 .5; PC1_cort_color];
 
 for i = 1:251
hold on
pIT(i).FaceColor=PC1_cort_color(i,:);
 end

print(['./GIF/PPEB_PC1.png'],'-dpng')
 
C1 = rescale([coeff{3}(1:250,1); -.04; .04],ylimits(1),ylimits(2))*-1;
C2 = rescale([IT; 5; 50],zlimits(1),zlimits(2));
 
for i =1:250
   [x{i},y{i}] = circle(C1(i),C2(i),2,length(BOUNDARY{i+1})); 
end

for i = 1:251
mean_expanded = mean(BOUNDARY{i}*1.33)-mean(BOUNDARY{i});

BOUNDARY_exploded{i} = BOUNDARY{i}+mean_expanded;
end


for j = 1:60
    shake = randINrange(-1,1,[1 3]);
    for i = 1:251
        pIT(i).Vertices = BOUNDARY{i}+shake;
    end
    ylim(ylimits)
zlim(zlimits)
xlim(xlimits)
pause(.1)
print(['./GIF/PPEB_S3_',num2str(j),'.png'],'-dpng')
end

r = 0:.1:1;
for j = 1:length(r)
for i = 1:251

bnew= find_point_on_line(BOUNDARY{i},BOUNDARY_exploded{i},r(j));
if i == 1
 pIT(i).FaceAlpha = 1-r(j);  
 pIT(i).EdgeAlpha = 1-r(j);  
end
hold on
pIT(i).Vertices = bnew;
material dull
end
ylim(ylimits)
zlim(zlimits)
xlim(xlimits)
pause(.1)
print(['./GIF/PPEB_S4_',num2str(j),'.png'],'-dpng')
end


r =0:.01:1;

for j = 1:length(r)
%r = .8;
for i = 1:250
  bnew= find_point_on_line(BOUNDARY_exploded{i+1},[ones(length(x{i}),1)*i x{i} y{i}],r(j));  
  pIT(i+1).FaceColor = find_point_on_line(IT_color(i+1,:),PC1_cort_color(i+1,:),r(j));  
pIT(i+1).Vertices = bnew;

end
ylim(ylimits)
zlim(zlimits)
xlim(xlimits)
%print(['TestParcel2Point',num2str(j),'.png'],'-dpng')
pause(.1)
print(['./GIF/PPEB_S5_',num2str(j),'.png'],'-dpng')
end

axis on

xaxis_pos = rescale(-.3:.1:.3,ylimits(1),ylimits(2));
xaxis_labels = -.04:.01:04;

yaxis_labels = 5:5:50;
yaxis_pos = rescale(yaxis_labels,zlimits(1),zlimits(2));

yticks(xaxis_pos)
yticklabels(xaxis_labels)

zticks(yaxis_pos)
zticklabels(yaxis_labels)

set(gca,'FontSize',24)
