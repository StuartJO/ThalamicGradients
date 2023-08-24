function PlotTractogramPCs()

% This script loads tractography data and associated principal component
% scores, and creates visualizations of tractograms based on these scores.
% It generates color-coded streamlines to represent tractography data,
% utilizing the first principal component (PC1) scores of both thalamic
% seed regions and cortical ROI coefficients. The resulting visualizations
% are saved as PNG images.
%
% The script performs the following steps:
% 1. Loads necessary data and calculates colormap mapping.
% 2. Creates a colormap image for visualization reference.
% 3. Creates a scatter plot of PC1 scores for thalamic seed regions and
%    cortical ROI coefficients.
% 4. Loops through each seed for example tracks and loads corresponding
%    track data and ROI assignments.
% 5. Calculates track lengths and filters out unwanted tracks.
% 6. Initializes arrays to store streamline coordinates and colors.
% 7. Loops through each seed again and processes track data, appending
%    streamline coordinates and colors.
% 8. Reshapes color data and creates a patch for streamlines.
% 9. Exports lateral and superior views of the tractogram visualizations.
%
% Output:
% The generated tractogram visualizations are saved as 'TractogramLateral.png'
% and 'TractogramSuperior.png' in the 'figure_outputs' directory.

% Load the necessary data
load('./data/processed/main_decomp.mat', 'pcs_thal', 'pcs_cort')

% Define colormap and mappings for tractogram visualization
numColors = 256;
ramp = linspace(-100, 100, numColors);
cform = makecform('lab2srgb');
a = repmat(ramp, [numColors 1]);             % -a on left
b = repmat(flipud(ramp'), [1 numColors]);    % -b on bottom
L = 50 * ones(numColors, numColors);         % A single L value.
Lab = cat(3, L, a, b);                       % A 2D image.
colormap2D = applycform(Lab, cform);

% Extract relevant data for visualization
PC_THAL = pcs_thal(:, 1);
PC_CORT = pcs_cort(:, 1);

% Rescale scores for colormap assignment
score_cmap_rescaled = rescale(zscore(PC_THAL), 1, 256);

% Rescale coefficients for colormap assignment
coeff_cmap_rescaled = rescale(zscore(PC_CORT), 1, 256);

% Create figure for colormap visualization
figure
imagesc(colormap2D);
set(gca, 'YDir', 'normal')
xlabel('Thalamic seed PC1 score')
ylabel('Cortical ROI PC1 coefficient')

% Configure tick labels for PC scores and coefficients
xticks_rescaled = rescale([-2:1 min(zscore(PC_THAL)) max(zscore(PC_THAL))], 1, 256);
yticks_rescaled = rescale([-2:1 min(zscore(PC_CORT)) max(zscore(PC_CORT))], 1, 256);
xticks(xticks_rescaled(1:end-2))
yticks(yticks_rescaled(1:end-2))
set(gca, 'XTickLabel', -2:1)
set(gca, 'YTickLabel', -2:1)
set(gca, 'FontSize', 20)
axis square

% Save colormap visualization as SVG and PNG
print('./figure_outputs/TractogramCmap.svg', '-dsvg', '-r300')
print('./figure_outputs/TractogramCmap.png', '-dpng', '-r300')

% Create a figure for tractogram visualization
figure('Position', [0 0 1920 1080])

% Configure lighting and view settings for 3D visualization
camlight(80, -10);     % Add light from upper-right direction
camlight(-80, -10);    % Add light from upper-left direction
view([90 0]);          % Set the initial view to a lateral view

% Configure axis settings for visualization
axis off;              % Turn off axis display
axis tight;            % Set the axis limits to fit the data
axis equal;            % Set equal scaling for all axis
axis vis3d;            % Set 3D visualization mode

% Create a waitbar to track progress
f = waitbar(0, 'Starting');

% Initialize cell array and counter for example tracks
MNIExampleTracks = cell(921, 3);
totalTractLengths = 0;

% Loop through all seeds for example tracks
for j = 1:921
    % Load example tracks and associated ROI assignments
    tracks = read_mrtrix_tracks(['./data/ExampleTracts/MNI/seed_', num2str(j), '.tck']);
    ROI = dlmread(['./data/ExampleTracts/conn/seed_', num2str(j), '_assignments']);
    tracks_data = tracks.data;
    Ntracks = length(tracks_data);
    Tracklen = zeros(length(tracks_data), 1);
    rmv = false(Ntracks, 1);
    
    % Calculate track lengths and identify tracks to remove
    for i = 1:Ntracks   
        Tracklen(i) = size(tracks_data{i}, 1);
        if ~(ROI(i) <= 250 && ROI(i) ~= 0)
            rmv(i) = true;        
        end
    end
    
    % Remove unwanted tracks and update data
    tracks_data(rmv) = [];
    ROI(rmv) = [];
    Tracklen(rmv) = [];
    MNIExampleTracks{j, 1} = tracks_data;
    MNIExampleTracks{j, 2} = ROI;
    MNIExampleTracks{j, 3} = Tracklen;
    
    % Update total tract lengths
    totalTractLengths = totalTractLengths + sum(Tracklen) + length(Tracklen);
    
    % Update waitbar progress
    waitbar(j / 921, f, ['Finished seed ', num2str(j), ' of ', num2str(921)]);
end
close(f)

% Initialize arrays to store streamlines and colors
STREAMLINES = zeros(totalTractLengths, 3);
STREAMLINES_COLOR = zeros(totalTractLengths, 3);

% Create another waitbar to track progress
f = waitbar(0, 'Starting');
IND = 0;

% Loop through all seeds again for streamline data
for j = 1:921
    tracks = MNIExampleTracks{j, 1};
    ROI = MNIExampleTracks{j, 2};
    
    % Loop through tracks associated with a seed
    for i = 1:length(tracks)
        tract = tracks{i};
        TractLength = length(tract) + IND;
        color2use = squeeze(colormap2D(round(score_cmap_rescaled(j)), round(coeff_cmap_rescaled(ROI(i))), :))';
        
        % Append streamline coordinates and colors
        STREAMLINES(IND+1:TractLength+1, :) = [[tract(:, 1)' NaN]', [tract(:, 2)' NaN]', [tract(:, 3)' NaN]'];
        STREAMLINES_COLOR(IND+1:TractLength+1, :) = repmat(color2use, size(tract, 1) + 1, 1);

        % Update index for next tract
        IND = TractLength + 1;
    end
    
    % Update waitbar progress
    waitbar(j / 921, f, ['Finished seed ', num2str(j), ' of ', num2str(921)]);
end
close(f)

% Reshape color data and create a patch for streamlines
cdata_all = reshape(STREAMLINES_COLOR, length(STREAMLINES_COLOR), 1, 3);    
streamlines = patch(STREAMLINES(:, 1), STREAMLINES(:, 2), STREAMLINES(:, 3), 0);
set(streamlines, 'CData', cdata_all, 'EdgeColor', 'interp', 'FaceColor', 'interp', 'Clipping', 'off') 

% Export lateral tractogram visualization
exportgraphics(gcf, './figure_outputs/TractogramLateral.png', 'Resolution', 300)

% Change the view to a superior view and export another visualization
view([90 90])
exportgraphics(gcf, './figure_outputs/TractogramSuperior.png', 'Resolution', 300)
