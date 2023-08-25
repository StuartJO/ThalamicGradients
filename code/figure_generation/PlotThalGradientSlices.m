function c = PlotThalGradientSlices(GRAD_DATA, GRAD_DATA_vox_coord, grad_cmap, cbar_title, thr, ImageOnly)

% Plot brain slices with thalamic gradient overlay.

% Inputs:
%   GRAD_DATA: Gradient data to overlay on the brain slices.
%   GRAD_DATA_vox_coord: Voxel coordinates of gradient data.
%   grad_cmap: Colormap for the gradient data.
%   cbar_title: Title for the plot/colorbar.
%   thr: Threshold for voxels to colour based on their distance to each seed.
%   ImageOnly: Flag to indicate whether to plot only the image without annotations or colorbar.

% Load necessary data
load('./data/preprocessed/MNI_Seed_voxelData.mat', 'braindata', 'thalmask')

% Check if ImageOnly input is provided, otherwise default to 0
if nargin < 6
    ImageOnly = 0;
end

% Find indices of voxels within the thalamic mask
inmask = find(thalmask == 1);

% Convert 1D indices to 3D voxel coordinates
[X, Y, Z] = ind2sub(size(thalmask), inmask);

% Compute distances between thalamic mask voxels and gradient data voxels
ThalMaskDists = pdist2([X Y Z], GRAD_DATA_vox_coord);

% Find the nearest gradient data voxel for each thalamic mask voxel
[ThalMask_nearest_seedmm, ThalMask_nearest_seed] = min(ThalMaskDists');

% Assign cluster IDs based on the nearest gradient data voxel
ThalVoxClust = changem(ThalMask_nearest_seed, GRAD_DATA, 1:length(GRAD_DATA_vox_coord));

% Create a thalamic mask with cluster IDs
thalmaskclust = thalmask;
thalmaskclust(thalmaskclust == 0) = NaN;
thalmaskclust(inmask) = ThalVoxClust;

% Mask out thalamic voxels based on distance threshold
thalmaskclust(inmask(ThalMask_nearest_seedmm > thr)) = NaN;

% Create a figure for plotting
figure

% Define custom colormap for slice annotations
slices_cmap = [254, 178, 76;...
    253, 141, 60;...
    252, 78, 42;...
    227, 26, 28;...
    189, 0, 38;...
    128, 0, 38]./255;

% Define the position of the brain slice axes
axes_slices = axes('Position', [.5 .57 .6 .5]);

% Display the brain slice
brain_back_slice = (squeeze(braindata(103, :, :)))';
imagesc(brain_back_slice);
set(gca, 'YDir', 'normal')
colormap(axes_slices, gray(256))

% Define slice positions
slices = [90 86 82 78 74 70];

% Plot lines for slice annotations if not in ImageOnly mode
if ~ImageOnly
    for i = 1:6
        hold on
        plot([1 size(brain_back_slice, 2)], [slices(i) slices(i)], 'Color', slices_cmap(i, :))
    end
end

% Set plot properties
axis image
axis off
axes_slices.XTick = [];
axes_slices.YTick = [];

% Define ranges for x and y slices
yslice_range = 87:127;
xslice_range = 60:120;

% Extract brain and thalamic slices
for i = 1:length(slices)
    slice = slices(i);
    brain_slices{i} = flipud(squeeze(braindata(xslice_range, yslice_range, slice))');
    thal_slices{i} = flipud(squeeze(thalmaskclust(xslice_range, yslice_range, slice))');
end

% Combine brain and thalamic slices
brain_back_slice = [brain_slices{1} brain_slices{2} brain_slices{3}; brain_slices{4} brain_slices{5} brain_slices{6}];
thal_slice = [thal_slices{1} thal_slices{2} thal_slices{3}; thal_slices{4} thal_slices{5} thal_slices{6}];

% Create axes for the background
axes_back = axes('Position', [0 0 1 .57]);

% Rescale brain and thalamic slices
brain_back_slice_rescaled = rescale(brain_back_slice, 1, 256);
thal_slice_rescaled = rescale(thal_slice, 257, 256 + size(grad_cmap, 1));

% Replace background pixels with thalamic gradient data
brain_back_slice_rescaled(~isnan(thal_slice)) = thal_slice_rescaled(~isnan(thal_slice));

% Create a combined colormap for background and gradient data
cmap = [gray(256); grad_cmap];

% Display the background image
imagesc(brain_back_slice_rescaled);
colormap(axes_back, cmap);
axes_back.XTick = [];
axes_back.YTick = [];
caxis(axes_back, [1 size(cmap, 1)])

% Add annotations and colorbar if not in ImageOnly mode
if ~ImageOnly
    XtextSpacing = size(brain_back_slice, 2) / 3;
    YtextSpacing = size(brain_back_slice, 1) / 2;
    x = [1 XtextSpacing (XtextSpacing * 2) 1 XtextSpacing (XtextSpacing * 2)] + 10;
    y = [1 1 1 YtextSpacing YtextSpacing YtextSpacing] + 20;

    for i = 1:6
        text(x(i), y(i), ['\bf', num2str(i)], 'Color', slices_cmap(i, :), 'FontSize', 28)
    end

    hold on
    plot([0 size(brain_back_slice, 2) + .5], [YtextSpacing YtextSpacing], 'k', 'LineWidth', 2.5)
    plot([XtextSpacing XtextSpacing], [0 size(brain_back_slice, 1) + .5], 'k', 'LineWidth', 2.5)
    plot([XtextSpacing XtextSpacing] * 2, [0 size(brain_back_slice, 1) + .5], 'k', 'LineWidth', 2.5)

    % Create axes for the gradient colorbar
    P = get(axes_back, 'pos');
    axes_front = axes('Position', P);
    set(axes_front, 'color', 'none')
    colormap(axes_front, grad_cmap);
    caxis(axes_front, [min(GRAD_DATA), max(GRAD_DATA)]);
    axes_front.XTick = [];
    axes_front.YTick = [];
    c = colorbar('Location', 'northoutside', 'FontSize', 18);

    % Adjust the position of the colorbar
    c.Position = [0.025, 0.638095238095239, 0.526785714285714, 0.066666666666666];

    % Add title annotation
    annotation(gcf, 'textbox',...
        [0.025, 0.754761904761905, 0.526785714285714, 0.238095238095242],...
        'VerticalAlignment', 'middle',...
        'String', cbar_title,...
        'HorizontalAlignment', 'center',...
        'FontSize', 18,...
        'FitBoxToText', 'off',...
        'EdgeColor', 'none');
end
