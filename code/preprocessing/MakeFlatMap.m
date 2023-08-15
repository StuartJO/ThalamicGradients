function MakeFlatMap()

load ('./data/preprocessed/AllenGeneDataset_19419.mat','structInfo')

% flat_map = flipud(fliplr(double(nrrdread('./data/ancillary/flatmap_dorsal.nrrd'))));
flat_map = rot90(double(nrrdread('./data/ancillary/flatmap_dorsal.nrrd')),2);

% FlatMapIds.xlsx was made by eyeballing Harris et al. Nature 2019, Fig 1b
% and comparing it to the initial boundaries drawn. I would use the values
% in flatmap_dorsal.nrrd to figure out what the regions but doing it by
% hand was faster. The IDs correspond to structure IDs not order IDs
FlatMapIds = readtable('./data/ancillary/FlatMapIds.xlsx');

% 'PTLp' isn't in this flat map but two of its subregions, 'VISa' and 
% 'VISrl' are. So lets just replace those
FlatMapIds_abbrev = [FlatMapIds.Abbreviation;'PTLp'];
replace2PTLpIND = ismember(FlatMapIds.Abbreviation,{'VISa','VISrl'});
replace2PTLp = FlatMapIds.Flatmap_id(replace2PTLpIND); 
Flatmap_id = [FlatMapIds.Flatmap_id;22];

Flatmap_id(replace2PTLpIND) = [];
FlatMapIds_abbrev(replace2PTLpIND) = [];

flat_map(flat_map==replace2PTLp(1)|flat_map==replace2PTLp(2))=22;

[~,I]=ismember(FlatMapIds_abbrev,structInfo.acronym);

% Set regions of flat map image that are empty space to be -1
I(1) = -1;
flat_map_matched = changem(flat_map,I,Flatmap_id);

flat_map_image_size = length(flat_map);

[X,Y] = meshgrid(1:flat_map_image_size);

% Downsample the image to make the boundaries look less irregular

[Xq,Yq] = meshgrid(2:2:flat_map_image_size);
Vq = interp2(X,Y,flat_map_matched,Xq,Yq,'nearest');

Vqsmooth = Vq;
for i = 2:length(Vq)-1
   for j = 2:length(Vq)-1
       nei = Vqsmooth(i-1:i+1,j-1:j+1);
       neiMode = mode(nei(:));
       if ~ismember(Vqsmooth(i,j),neiMode)
           Vqsmooth(i,j) = neiMode(1);
       end
   end
end

% Convert to a triangular mesh to make pretty boundaries

[F,V]=mesh2tri(Xq,Yq,Vqsmooth,'f');

vertex_id = double(V(:,3));
vertices = double([V(:,1:2) zeros(size(V,1),1)]);
faces = double(F);

[Mouseflatmap_boundary,Mouseflatmap_boundary_id] = findROIboundaries(vertices,faces,vertex_id,'midpoint');

% Again, eyeballing Harris et al. Nature 2019, Fig 1f to get this. It
% sucked

MouseCorticalDivisions = [2     6     6     6     5     5     5     6     1     1     1     1     1     1     1     4 ...
    2     4     4     4     2     5     5     5     6     6     2     2     2     5     5     6 ...
    1     1     3     3     3     5];

save('./data/ancillary/Mouseflatmap.mat','Mouseflatmap_boundary','Mouseflatmap_boundary_id','MouseCorticalDivisions')