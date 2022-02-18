function [ObsClust,VarClust] = PlotClustergram(DATA,NObsClust,NVarClust,cmap)

Nobs = size(DATA,1);
Nvar = size(DATA,2);

MatPos = [0.15 0.01 0.6604 0.7495];

matAx = axes('Position',MatPos);

YdenAx_xpos = .05;

YdenAx = axes('Position',[YdenAx_xpos MatPos(2) MatPos(1)-YdenAx_xpos MatPos(4)]);

XdenAx_ypos = MatPos(2)+MatPos(4);

XdenAx = axes('Position',[MatPos(1) XdenAx_ypos MatPos(3) .1]);

D_y = squareform(pdist(DATA));

Z_y = linkage(DATA,'ward');

order_y = optimalleaforder(Z_y,D_y);

D_x = squareform(pdist(DATA'));

Z_x = linkage(DATA','ward');

order_x = optimalleaforder(Z_x,D_x);

data = DATA(flip(order_y),order_x);
%data = D(flip(order_x),flip(order_x));
axes(matAx)
imagesc(matAx,data)
axis off

SeedClust = NObsClust;
ObsClust = cluster(Z_y,'maxclust',SeedClust);
%cutoff = median([Z_y(end-SeedClust+1,3) Z_y(end-SeedClust+2,3)]);
axes(YdenAx)
%h_obs = dendrogram(Z_y,0,'reorder',order_y,'Orientation','Left','ColorThreshold',cutoff);
h_obs = dendrogram(Z_y,0,'reorder',order_y,'Orientation','Left');
set(h_obs,'Color',[0 0 0])
ylim([.5 Nobs+.5])
xlimits = xlim;
if NObsClust ~= 1
for i = 1:NObsClust
clust_x = find(ObsClust(order_y)==i);
clust_start = min(clust_x)-.5;
clust_range =  max(clust_x)-clust_start+.5;
rectangle('Position',[0,clust_start,xlimits(2),clust_range],'FaceColor',cmap(i,:),'EdgeColor','None')
end
chi=get(gca, 'Children');
set(gca, 'Children',flipud(chi));
end
axis off

yyaxis left
yticks([])
ylabel(YdenAx,'Thalamus seeds','FontWeight','bold')
YdenAx.YLabel.Visible = 'on';

NodeClust = NVarClust;
VarClust = cluster(Z_x,'maxclust',NodeClust);
%cutoff = median([Z_y(end-NodeClust+1,3) Z_y(end-NodeClust+2,3)]);
axes(XdenAx)
%h_vars = dendrogram(Z_x,0,'reorder',order_x,'ColorThreshold',cutoff);
h_vars = dendrogram(Z_x,0,'reorder',order_x);
set(h_vars,'Color',[0 0 0])
xlim([.5 Nvar+.5])
ylimits = ylim;
if NVarClust ~= 1
for i = 1:NVarClust
clust_x = find(VarClust(order_x)==i);
clust_start = min(clust_x)-.5;
clust_range =  max(clust_x)-clust_start+.5;
rectangle('Position',[clust_start,0,clust_range,ylimits(2)],'FaceColor',cmap(i,:),'EdgeColor','None');
end
chi=get(gca, 'Children');
set(gca, 'Children',flipud(chi));
end
axis off
title('Cortical regions')
axes(matAx)



