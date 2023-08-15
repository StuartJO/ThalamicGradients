function GetCorticalSpinTestPerms()
% Computes spin test permutations for cortical parcellation
%
%   GetCorticalSpinTestPerms()
%
%   This function computes spin test permutations for a cortical parcellation
%   based on spherical surface data. The permutations are used to assess the
%   statistical significance of observed patterns or relationships in the data.

% Define the number of iterations for spin test permutations
Niters = 10000;

% Load the left hemisphere spherical surface
lh_sphere = read_surface('lh.sphere');

% Load randomly selected vertices on the fsaverage surface
load('fsaverage_surface_data.mat', 'lh_rand500')

% Initialize an array to store the centroids of cortical regions
centroid = zeros(250, 3);

% Calculate the centroid coordinates for each cortical region
for i = 1:250
    centroid(i, :) = mean(lh_sphere.vertices(lh_rand500 == i, :));
end

% Perform spin test permutations to create permuted parcellations
perm_id = rotate_parcellation_single(centroid, Niters);

% Save the generated spin test permutations
save('./data/preprocessed/CorticalSpinTestPerms.mat', 'perm_id')
