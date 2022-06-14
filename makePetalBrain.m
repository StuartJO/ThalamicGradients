figure('Position',[0 0 1920 1080])

load('IT.mat')
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

% camlight(80,-10);
% camlight(-80,-10);
view([-90 0])
axis off
axis image

ylimits = ylim;
zlimits = zlim;
xlimits = xlim;

for i = 1:251
hold on
pIT(i).FaceColor=[.5 .5 .5];
end

print(['./NewGIFs/PPEB_S1.png'],'-dpng')
view([90 0])
print(['./NewGIFs/PPEB_S1_medial.png'],'-dpng')

 PC1_color = turbo(256);
 PC1_cort_color = value2color(zscore(coeff{3}(1:250,1)),PC1_color);
 
 PC1_cort_color = [.5 .5 .5; PC1_cort_color];
 
 for i = 1:251
hold on
pIT(i).FaceColor=PC1_cort_color(i,:);
 end
 
view([-90 0])
print(['./NewGIFs/PPEB_PC1.png'],'-dpng')
view([90 0])
print(['./NewGIFs/PPEB_PC1_medial.png'],'-dpng')

for i = 1:251
hold on
pIT(i).FaceColor=IT_color(i,:);
end
view([-90 0])
print(['./NewGIFs/PPEB_S2.png'],'-dpng')

 
C1 = rescale([zscore(coeff{3}(1:250,1)); -2.5; 1.5],ylimits(1),ylimits(2))*-1;
C2 = rescale([IT; 5; 50],zlimits(1),zlimits(2));
 
for i =1:250
   [x{i},y{i}] = circle(C1(i),C2(i),2,length(BOUNDARY{i+1})); 
end

for i = 1:251
   meanYpos(i) = mean(BOUNDARY{i}(:,2));
end

[~,Order] = sort(meanYpos,'descend');

%AnchorPoint = meanYpos+10;

R = zeros(251,60);

r =linspace(0,1,30);

for i = 1:251
    StartTime = round(30*(find(Order==i)/251))+1;
    R(i,StartTime:StartTime+29) = r;
    R(i,StartTime+30:end) = 1;
    if i == 1
      medWall = zeros(60,1);
      medWall(1:StartTime) = linspace(1,0,StartTime);
    end
end


for j = 1:size(R,2)
%r = .8;
for i = 1:251
    if i == 1
     pIT(i).FaceAlpha = medWall(j);  
     pIT(i).EdgeAlpha = medWall(j);  
    else
    
    
  bnew= find_point_on_line(BOUNDARY{i},[ones(length(x{i-1}),1)*(Order(i)) x{i-1} y{i-1}],R(i,j));  
  pIT(i).FaceColor = find_point_on_line(IT_color(i,:),PC1_cort_color(i,:),R(i,j));  
pIT(i).Vertices = bnew;
    end
end
ylim(ylimits)
zlim(zlimits)
xlim(xlimits)
%print(['TestParcel2Point',num2str(j),'.png'],'-dpng')
pause(.1)
print(['./NewGIFs/PPEB_S3_',num2str(j),'.png'],'-dpng')
end

axis on

xaxis_pos = rescale(-2.5:.5:1.5,ylimits(1),ylimits(2));
xaxis_labels = 1.5:-.5:-2.5;

yaxis_labels = 5:5:50;
yaxis_pos = rescale(yaxis_labels,zlimits(1),zlimits(2));

yticks(xaxis_pos)
yticklabels(xaxis_labels)

zticks(yaxis_pos)
zticklabels(yaxis_labels)

set(gca,'FontSize',24)
print(['./NewGIFs/PPEB_END.png'],'-dpng')