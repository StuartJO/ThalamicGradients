function MouseThalROICoords = GetMouseThalROIcentroid(Simple)
% Computes centroid coordinates of thalamic regions
%
%   MouseThalROICoords = GetMouseThalROIcentroid(Simple)
%
%   Inputs:
%   - Simple: A logical flag indicating whether to compute simple centroids (default: false).
%             If Simple is true, simple centroids (center of gravity) are computed for each
%             thalamic region. If false, the point within the region that minimizes the mean
%             distance to all other points in the region is used as the centroid.
%
%   Output:
%   - MouseThalROICoords: A 35x3 numeric array containing the centroid coordinates (x, y, z)
%                         for each thalamic region.
%
%   This function calculates the centroid coordinates for each thalamic region in the Allen Mouse
%   Brain Atlas based on the provided thalamic parcellation data. The centroids can be computed
%   either as simple centroids (center of gravity) or as points that minimize mean distances.
%   The computed coordinates are returned in the 'MouseThalROICoords' array.

% Set default value for Simple if not provided
if nargin < 1
    Simple = 0;
end

% Load thalamic parcellation data
load('./data/ancillary/MouseOhParc.mat', 'MouseOhParc')

% Define thalamic regions of interest
ThalRegions = 88:122;

% Prepare thalamic parcellation data
MouseThalOnly = MouseOhParc;
MouseThalOnly(~ismember(MouseOhParc, ThalRegions)) = NaN;
MouseThalOnly(MouseThalOnly == 0) = NaN;
MouseThalOnly(1:228, :, :) = NaN;

% Create thalamic parcellation with region labels
MouseThalROI = changem(MouseThalOnly, 1:35, ThalRegions);

% Initialize array to store centroid coordinates
MouseThalROICoords = zeros(35, 3);

% Loop through thalamic regions
for i = 1:35
    % Find indices of voxels within the region
    IND = find(MouseThalROI == i);
    [mX, mY, mZ] = ind2sub(size(MouseThalROI), IND);
    
    % Compute centroid based on chosen method
    if Simple
        MouseThalROICoords(i, :) = mean([mX, mY, mZ]); % Simple centroid (center of gravity)
    else
        d = squareform(pdist([mX, mY, mZ]));
        [~, I] = min(mean(d));
        MouseThalROICoords(i, :) = [mX(I), mY(I), mZ(I)]; % Centroid minimizing mean distances
    end
    
    disp(['Got x, y, z coords for thalamic region ', num2str(i)])
end

% Save the computed coordinates based on chosen method
if Simple == 0
    save('./data/ancillary/MouseThalROICoords.mat', 'MouseThalROICoords')
else
    save('./data/ancillary/MouseThalROICoords_COG.mat', 'MouseThalROICoords')
end