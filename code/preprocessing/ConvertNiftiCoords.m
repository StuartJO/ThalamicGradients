function outcoords = ConvertNiftiCoords(nifti_info, coords, inCoordType, VoxStartIND)

% Convert coordinates between voxel and millimeter space using the provided 
% transformation matrix.

% Inputs:
%   nifti_info: NIFTI header information containing the transformation 
%   matrix. This can be obtain with the "niftiinfo" command
%   inCoordType: Coordinates to be converted. Should be a matrix with each row representing a set of coordinates.
%   intype: Input coordinate space type ('mm' for millimeters or 'vox' for voxels).
%   VoxStartIND: Flag indicating if voxel indexing starts from 1 (default) 
%   or 0. If working in MATLAB you'll likely want it to begin at 1, but if
%   importing/exporting coordinates from external neuroimaging software 
%   (like FSL or MRtrix), you'll want to star at 0
%
% Outputs:
%   outcoords: Coordinates in the opposite space to the input (i.e., if the
%   input were in mm, then the output iwll be voxels
%
% Set a default value for VoxStartIND if not provided
if nargin < 4
    VoxStartIND = 1;
end

% Extract the transformation matrix from nifti_info
T = nifti_info.Transform.T;

% Compute the voxel-to-millimeter and millimeter-to-voxel transformation matrices
v2m = T';
voxscale = diag(v2m(1:3, 1:3));
if VoxStartIND == 1
    v2m(1:3, 4) = v2m(1:3, 4) - voxscale;
end
m2v = inv(v2m);

% Initialize the output coordinates matrix
outcoords = zeros(size(coords));

% Switch between input coordinate types
switch inCoordType
    case 'mm'
        mm = coords;
        % Convert millimeter coordinates to voxel coordinates
        for i = 1:size(mm, 1)
            outcoords(i, 1:3) = mm(i, :) * m2v(1:3, 1:3) + m2v(1:3, 4)';
        end    
    case 'vox'
        vox = coords;
        % Convert voxel coordinates to millimeter coordinates
        for i = 1:size(vox, 1)
            outcoords(i, :) = (vox(i, 1:3) - m2v(1:3, 4)') / m2v(1:3, 1:3);
        end        
end