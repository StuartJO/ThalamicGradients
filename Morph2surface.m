figure('Position',[0 0 1920 1080])
addpath F:\Documents\GitHub\plotSurfaceROIBoundary
axes('Position',[0.0177 0.0241 0.9698 0.9657])

load('IT.mat')
IT_colormap = cbrewer('div','RdBu',256);

IT_color_ind = ceil(rescale(IT,1,size(IT_colormap,1)));
IT_color = IT_colormap(IT_color_ind,:);
IT_color = [.5 .5 .5; IT_color];

 PC1_color = turbo(256);
 PC1_cort_color = value2color(zscore(coeff{3}(1:250,1)),PC1_color);
 PC1_cort_color = [.5 .5 .5; PC1_cort_color];

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

% camlight(80,-10);
% camlight(-80,-10);
view([-90 0])
axis off
axis image

ylimits = ylim;
zlimits = zlim;
xlimits = xlim;

print(['./NewGIFs/Morph2surface_START.png'],'-dpng')

plotmiddle = mean(ylimits);

% z_range = linspace(120,40,251);
% y_range = [plotmiddle-20 plotmiddle+20];

z_range = linspace(120,-120,251);
y_range = [plotmiddle-25 plotmiddle+25];
%y_rec = [y_range(1) y_range(2)-y_range(1)];

y_rec1 = [y_range(1) (y_range(2)-y_range(1))/2];
y_rec2 = [y_range(1)+(y_range(2)-y_range(1)) -(y_range(2)-y_range(1))/2];


for i = 1:250
    b = BOUNDARY{i+1};
    npoints = size(b,1);
    z_rec = [z_range(i+1) z_range(i)-z_range(i+1)];
%     [x_,y_] = rectangleNpoints(y_rec,z_rec,npoints,1);
%     rec_coords{i} = [ones(npoints,1)*-50 x_ y_];
    [x_,y_] = rectangleNpoints(y_rec1,z_rec,npoints,1);
    rec_coords_med{i} = [ones(npoints,1)*-70 x_ y_];
    [x_,y_] = rectangleNpoints(y_rec2,z_rec,npoints,1);
    rec_coords_lat{i} = [ones(npoints,1)*-70 x_ y_];
end

hold on

for i = 1:250
hold on
pMed(i)=fill3(rec_coords_med{i}(:,1),rec_coords_med{i}(:,2),rec_coords_med{i}(:,3),PC1_cort_color(i+1,:),'FaceLighting','none','Clipping','off','LineStyle','none');
pLat(i)=fill3(rec_coords_lat{i}(:,1),rec_coords_lat{i}(:,2),rec_coords_lat{i}(:,3),PC1_cort_color(i+1,:),'FaceLighting','none','Clipping','off','LineStyle','none');

material dull

end

[~,Order] = sort(zscore(coeff{3}(1:250,1)),'descend');

frames = 45;
N2target = 10;
R = zeros(250,frames);
r =linspace(0,1,N2target);

last_time = frames-N2target;

for i = 1:250
    StartTime = round(last_time*(find(Order==i)/250))+1;
    R(i,StartTime:StartTime+N2target-1) = r;
    R(i,StartTime+N2target:end) = 1;
end

for j = 1:size(R,2)
%r = .8;
for i = 1:250
       
    bnew= find_point_on_line(rec_coords_med{i},BOUNDARY{i+1},R(i,j));    
    pMed(i).Vertices = bnew;
    
    bnew= find_point_on_line(rec_coords_lat{i},BOUNDARY2{i+1},R(i,j));    
    pLat(i).Vertices = bnew;
    
    if R(i,j) == 1
        pMed(i).FaceAlpha = 0;
        pIT(i+1).FaceColor = PC1_cort_color(i+1,:);
        
        pLat(i).FaceAlpha = 0;
        pIT2(i+1).FaceColor = PC1_cort_color(i+1,:);
    else
        pMed(i).FaceAlpha = 1;
        pIT(i+1).FaceColor = [.5 .5 .5];
        
        pLat(i).FaceAlpha = 1;
        pIT2(i+1).FaceColor = [.5 .5 .5];
    end
end
ylim(ylimits)
zlim(zlimits)
xlim(xlimits)
%print(['TestParcel2Point',num2str(j),'.png'],'-dpng')
pause(.1)
print(['./NewGIFs/Morph2surface2_',num2str(j),'.png'],'-dpng')
end

for i = 1:251
    pIT(i).FaceColor = IT_color(i,:);
    pIT2(i).FaceColor = IT_color(i,:);
end
print(['./NewGIFs/Morph2surface_IT.png'],'-dpng')