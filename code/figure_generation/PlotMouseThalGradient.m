function PlotMouseThalGradient(ThalData, BackGroundData, grad_cmap, cbar_title)
% PlotMouseThalGradient Plots thalamic gradient data on brain slices
%
%   PlotMouseThalGradient(ThalData, BackGroundData, grad_cmap, TITLE)
%
%   Inputs:
%   - ThalData: Thalamic gradient data to be plotted.
%   - BackGroundData: Background brain data for visualization.
%   - grad_cmap: Colormap used for visualization of thalamic gradient data.
%   - cbar_title: Title for the plot.
%
%   This function generates a visualization of thalamic gradient data overlaid on
%   brain slices. It creates a figure with labeled slices and colormap.

% Create a new figure
figure

% Define slices colormap for labeling
slices_cmap = [254,178,76;...
    253,141,60;...
    252,78,42;...
    227,26,28;...
    189,0,38;...
    128,0,38]./255;

% Create axes for slices
axes_slices = axes('Position', [0.5, 0.57, 0.5, 0.43]);
brain_back_slice = (squeeze(BackGroundData(220, :, :))');
h_brain_back = imagesc(brain_back_slice);
set(gca, 'YDir', 'normal')
colormap(axes_slices, gray(256))

% Define slices and draw lines for labeling
slices = round(linspace(140, 180, 6));
for i = 1:6
    hold on
    plot([1, size(brain_back_slice, 2)], [slices(i), slices(i)], 'Color', slices_cmap(i, :))
end
axes_slices.XTick = [];
axes_slices.YTick = [];

% Define ranges for slice extraction
yslice_range = 160:310;
xslice_range = 175:385;

% Extract slices from background and thalamic data
for i = 1:length(slices)
    slice = slices(i);
    brain_slices{i} = flipud(squeeze(BackGroundData(xslice_range, yslice_range, slice))');
    thal_slices{i} = flipud(squeeze(ThalData(xslice_range, yslice_range, slice))');
end

% Combine slices into complete images
brain_back_slice = [brain_slices{1}, brain_slices{2}, brain_slices{3}; brain_slices{4}, brain_slices{5}, brain_slices{6}];
thal_slice = [thal_slices{1}, thal_slices{2}, thal_slices{3}; thal_slices{4}, thal_slices{5}, thal_slices{6}];

% Create axes for background and overlay
axes_back = axes('Position', [0, 0, 1, 0.57]);

% Rescale background and thalamic data
brain_back_slice_rescaled = rescale(brain_back_slice, 1, 256);
thal_slice_rescaled = rescale(thal_slice, 257, 256 + size(grad_cmap, 1));

% Replace background pixels with thalamic data
brain_back_slice_rescaled(~isnan(thal_slice)) = thal_slice_rescaled(~isnan(thal_slice));

% Create colormap with transparency
gray_cmap = make_alpha_rgb(gray(256), 0.2);
cmap = [gray_cmap; grad_cmap];

% Display background and overlay data
h_brain_back = imagesc(brain_back_slice_rescaled);
colormap(axes_back, cmap);

axes_back.XTick = [];
axes_back.YTick = [];
caxis(axes_back, [1, size(cmap, 1)])

% Define text positions for slice numbers
XtextSpacing = size(brain_back_slice, 2) / 3;
YtextSpacing = size(brain_back_slice, 1) / 2;
x = [1, XtextSpacing, (XtextSpacing * 2), 1, XtextSpacing, (XtextSpacing * 2)] + 10;
y = [1, 1, 1, YtextSpacing, YtextSpacing, YtextSpacing] + 20;

% Place slice numbers on the plot
for i = 1:6
    text(x(i), y(i), ['\bf', num2str(i)], 'Color', slices_cmap(i, :), 'FontSize', 28)
end
hold on
plot([0, size(brain_back_slice, 2) + 0.5], [YtextSpacing, YtextSpacing], 'k', 'LineWidth', 2) 
plot([XtextSpacing, XtextSpacing], [0, size(brain_back_slice, 1) + 0.5], 'k', 'LineWidth', 2) 
plot([XtextSpacing, XtextSpacing] * 2, [0, size(brain_back_slice, 1) + 0.5], 'k', 'LineWidth', 2) 

% Create a new set of axes for the colormap
P = get(axes_back, 'pos');
axes_front = axes('Position', P);
set(axes_front, 'color', 'none')
colormap(axes_front, grad_cmap);

% Set color axis range
caxis(axes_front, [nanmin(ThalData(:)), nanmax(ThalData(:))]);

axes_front.XTick = [];
axes_front.YTick = [];

% Add colorbar at the top of the plot
c = colorbar('Location', 'northoutside', 'FontSize', 18);
c.Position = [0.025, 0.638095238095239, 0.5 - 0.05, 0.066666666666666];

% Add colorbar title annotation
annotation(gcf, 'textbox',...
    [0, 0.754761904761905, 0.5, 0.238095238095242],...
    'VerticalAlignment', 'middle',...
    'String', cbar_title,...
    'HorizontalAlignment', 'center',...
    'FontSize', 18,...
    'FitBoxToText', 'off',...
    'EdgeColor', 'none');
