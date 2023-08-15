function MappedColors = MapData2Colors(data, cmap, limits)
% MapData2Colors: A function to map data values to colors using a colormap
%
% Inputs:
%   data  - Data values to be mapped to colors
%   cmap  - Colormap to use for mapping (e.g., obtained from 'brewermap' function)
%   limits (optional) - Limits for mapping data values to colors (default: min and max of data)
%
% Outputs:
%   MappedColors - RGB colors corresponding to the mapped data values
%
% Usage example:
%   mapped_colors = MapData2Colors(data, cmap);
%
%   This function maps data values to colors using the provided colormap and
%   ensures that the data values are scaled within the colormap limits. NaN
%   values are handled by mapping them to the first color in the colormap.

% Check if limits were provided, else set limits to the min and max of the data
    if nargin < 3
        limits = [min(data), max(data)];
    end

    % Prepare data and limits for mapping
    if size(data, 2) == 1
        data_ = [data; limits'];
    else
        data_ = [data, limits];
    end

    % Calculate color indices by scaling data to fit within the colormap
    color_ind = ceil(rescale(data_, 1, size(cmap, 1)));
    
    % Handle NaN values by assigning them to the first color in the colormap
    color_ind(isnan(color_ind)) = 1;

    % Use the color indices to retrieve corresponding colors from the colormap
    MappedColors = cmap(color_ind(1:end - 2), :);
end
