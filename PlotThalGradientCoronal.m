%function PlotThalGradientCoronal(GRAD_DATA,voxel_seeds,grad_cmap,TITLE,thr)

%grad_cmap = parula(256);
GRAD_DATA=pc_thal;
voxel_seeds = seed_voxel_coords;
grad_cmap = turbo(256);
TITLE=['Thalamic seed PC1 score'];
thr = 2.1;
%PlotThalGradient3(pc_thal,seed_voxel_coords,turbo(256),['Thalamic seed PC',num2str(grad),' score'],2.1)

load('Seed_voxelData.mat')

inmask = find(thalmask==1);
[X,Y,Z] = ind2sub(size(thalmask),inmask);
ThalMaskDists = pdist2([X Y Z],voxel_seeds);
[ThalMask_nearest_seedmm,ThalMask_nearest_seed] = min(ThalMaskDists');

ThalVoxClust = changem(ThalMask_nearest_seed,GRAD_DATA,1:length(voxel_seeds));

thalmaskclust = thalmask;
thalmaskclust(thalmaskclust==0) = NaN;
thalmaskclust(inmask) = ThalVoxClust;

thalmaskclust(inmask(ThalMask_nearest_seedmm>thr)) = NaN;

run = 1;

if run

figure

slices_cmap = [254,178,76;...
253,141,60;...
252,78,42;...
227,26,28;...
189,0,38;...
128,0,38]./255;
axes_slices = axes('Position',[.5 .57 .6 .5]);
brain_back_slice = (squeeze(brain(103,:,:))');
h_brain_back = imagesc(brain_back_slice);
set(gca,'YDir','normal')
colormap(axes_slices,gray(256))
slices = [90 84 78 72 66 60]+35;
for i = 1:6
    hold on
plot([slices(i) slices(i)],[1 size(brain_back_slice,1)],'Color',slices_cmap(i,:))
end
axis image
axes_slices.XTick = [];
axes_slices.YTick = [];

xslice_range = 60:120;
zslice_range = 65:95;
for i = 1:length(slices)
slice = slices(i);
brain_slices{i} = flipud(squeeze(brain(xslice_range,slice,zslice_range))');
thal_slices{i} = flipud(squeeze(thalmaskclust(xslice_range,slice,zslice_range))');
end

brain_back_slice = [brain_slices{1} brain_slices{2} brain_slices{3}; brain_slices{4} brain_slices{5} brain_slices{6}];
thal_slice = [thal_slices{1} thal_slices{2} thal_slices{3}; thal_slices{4} thal_slices{5} thal_slices{6}];

axes_back = axes('Position',[0 0 1 1]);
h_brain_back = imagesc(brain_back_slice);
colormap(axes_back,gray(256));
axis image

P = get(axes_back,'pos');
axes_front = axes('Position',P);
h = imagesc(thal_slice);
set(h, 'AlphaData', ~isnan(thal_slice))
set(axes_front,'color','none')
axis image
axes_front.XTick = [];
axes_front.YTick = [];
axes_back.XTick = [];
axes_back.YTick = [];
caxis([min(GRAD_DATA) max(GRAD_DATA)]);
colormap(axes_front,grad_cmap);
axes_back.Position = [0 -.2 1 1];
axes_front.Position = [0 -.2 1 1];

XtextSpacing = size(brain_back_slice,2)/3;
YtextSpacing = size(brain_back_slice,1)/2;

x = [1 XtextSpacing (XtextSpacing*2) 1 XtextSpacing (XtextSpacing*2)]+10;
y = [1 1 1 YtextSpacing YtextSpacing YtextSpacing]+10;

for i = 1:6
text(x(i),y(i),['\bf',num2str(i)],'Color',slices_cmap(i,:),'FontSize',28)
end
hold on
plot([0 size(brain_back_slice,2)+.5],[YtextSpacing YtextSpacing],'k','LineWidth',2) 

plot([XtextSpacing XtextSpacing],[0 size(brain_back_slice,1)+.5],'k','LineWidth',2) 
plot([XtextSpacing XtextSpacing]*2,[0 size(brain_back_slice,1)+.5],'k','LineWidth',2) 

c = colorbar('Location','northoutside','FontSize',18);

c.Position = [0.025,0.638095238095239,0.526785714285714,0.066666666666666];


annotation(gcf,'textbox',...
    [0,0.754761904761905,0.557142857142857,0.238095238095242],...
    'VerticalAlignment','middle',...
    'String',TITLE,...
    'HorizontalAlignment','center',...
    'FontSize',18,...
    'FitBoxToText','off',...
    'EdgeColor','none');

end