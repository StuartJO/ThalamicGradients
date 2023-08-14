function PlotMouseThalGradientAlt(ThalData,BackGroundData,grad_cmap,TITLE)

figure

slices_cmap = [254,178,76;...
253,141,60;...
252,78,42;...
227,26,28;...
189,0,38;...
128,0,38]./255;
axes_slices = axes('Position',[.5 .57 .5 .43]);
brain_back_slice = (squeeze(BackGroundData(220,:,:))');
h_brain_back = imagesc(brain_back_slice);
set(gca,'YDir','normal')
colormap(axes_slices,gray(256))
slices = round(linspace(140,180,6));
for i = 1:6
    hold on
plot([1 size(brain_back_slice,2)],[slices(i) slices(i)],'Color',slices_cmap(i,:))
end
%axis image
axes_slices.XTick = [];
axes_slices.YTick = [];

yslice_range = 160:310;
xslice_range = 175:385;

%yslice_range = 180:320;
%xslice_range = 195:405;

for i = 1:length(slices)
slice = slices(i);
brain_slices{i} = flipud(squeeze(BackGroundData(xslice_range,yslice_range,slice))');
thal_slices{i} = flipud(squeeze(ThalData(xslice_range,yslice_range,slice))');
end


brain_back_slice = [brain_slices{1} brain_slices{2} brain_slices{3}; brain_slices{4} brain_slices{5} brain_slices{6}];
thal_slice = [thal_slices{1} thal_slices{2} thal_slices{3}; thal_slices{4} thal_slices{5} thal_slices{6}];

axes_back = axes('Position',[0 0 1 .57]);

brain_back_slice_rescaled = rescale(brain_back_slice,1,256);

thal_slice_rescaled = rescale(thal_slice,257,256+size(grad_cmap,1));

brain_back_slice_rescaled(~isnan(thal_slice)) = thal_slice_rescaled(~isnan(thal_slice));

gray_cmap = make_alpha_rgb(gray(256),.2);

cmap = [gray_cmap; grad_cmap];

h_brain_back = imagesc(brain_back_slice_rescaled);
colormap(axes_back,cmap);

axes_back.XTick = [];
axes_back.YTick = [];

caxis(axes_back,[1 size(cmap,1)])

XtextSpacing = size(brain_back_slice,2)/3;
YtextSpacing = size(brain_back_slice,1)/2;

x = [1 XtextSpacing (XtextSpacing*2) 1 XtextSpacing (XtextSpacing*2)]+10;
y = [1 1 1 YtextSpacing YtextSpacing YtextSpacing]+20;

for i = 1:6
text(x(i),y(i),['\bf',num2str(i)],'Color',slices_cmap(i,:),'FontSize',28)
end
hold on
plot([0 size(brain_back_slice,2)+.5],[YtextSpacing YtextSpacing],'k','LineWidth',2) 

plot([XtextSpacing XtextSpacing],[0 size(brain_back_slice,1)+.5],'k','LineWidth',2) 
plot([XtextSpacing XtextSpacing]*2,[0 size(brain_back_slice,1)+.5],'k','LineWidth',2) 

P = get(axes_back,'pos');
axes_front = axes('Position',P);
set(axes_front,'color','none')
colormap(axes_front,grad_cmap);

caxis(axes_front,[nanmin(ThalData(:)) nanmax(ThalData(:))]);

axes_front.XTick = [];
axes_front.YTick = [];

c = colorbar('Location','northoutside','FontSize',18);

c.Position = [0.025,0.638095238095239,.5-0.05,0.066666666666666];

annotation(gcf,'textbox',...
    [0,0.754761904761905,0.5,0.238095238095242],...
    'VerticalAlignment','middle',...
    'String',TITLE,...
    'HorizontalAlignment','center',...
    'FontSize',18,...
    'FitBoxToText','off',...
    'EdgeColor','none');