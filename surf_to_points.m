

vertices = surface.vertices;
faces = surface.faces;
parc = lh_rand500;

parc_ids = unique(parc);

[~,~,~,BOUNDARY_VERTS] = findROIboundaries(vertices,faces,parc,'edge_vertices');

[BOUNDARY,~,~,~] = findROIboundaries(vertices,faces,parc,'midpoint');

for i =1:251
   [x{i},y{i}] = circle(C1,C2,2,length(BOUNDARY{i})); 
end

for i = 1:length(parc_ids)
    ROI_verts_ind = find(parc==parc_ids(i));
    ROI_verts = vertices(ROI_verts_ind,:);
faces_IND = find(sum(ismember(faces,ROI_verts_ind),2) > 0);

nROI_verts = size(ROI_verts,1);

ROIfaces = double(faces(faces_IND,:));

ROIfaces = ROIfaces.*ismember(ROIfaces,ROI_verts_ind);

%ROIfaces(sum(ROIfaces==0,2)>0,:) = [];

ROI_verts_new = changem(ROI_verts_ind,1:nROI_verts,ROI_verts_ind);

ROIfaces = changem(ROIfaces,1:nROI_verts,ROI_verts_ind);

ROI_edgeList = [[ROIfaces(:,1); ROIfaces(:,1); ROIfaces(:,2); ROIfaces(:,2); ROIfaces(:,3); ROIfaces(:,3)], ...
    [ROIfaces(:,3); ROIfaces(:,2); ROIfaces(:,1); ROIfaces(:,3); ROIfaces(:,2); ROIfaces(:,1)]];

ROI_edgeList(sum(ROI_edgeList==0,2)>0,:)=[];

adj = sparse(ROI_edgeList(:, 1), ROI_edgeList(:, 2), 1, nROI_verts, nROI_verts);
g = graph(adj);

[ROIboundary_verts,ROIboundary_verts_order] = ismember(ROI_verts_ind,BOUNDARY_VERTS{i}(1:end));
ROIboundary_verts_ind = find(ROIboundary_verts);

ROIboundary_verts_order = ROIboundary_verts_order(ROIboundary_verts_ind);

d = distances(g);

[~,middleVERT] = min(mean(d));

[closest_Boundary_vert_dist,closest_Boundary_vert_ind] = min(d(ROIboundary_verts_ind,:));

closest_Boundary_vert = ROIboundary_verts_ind(closest_Boundary_vert_ind);

[~,closestBOUNDARYpoint] = min(pdist2(ROI_verts(ROIboundary_verts_ind,:),BOUNDARY{i}),[],2);

BoundaryROIsend_coords = [ones(length(closestBOUNDARYpoint),1) x{i}(closestBOUNDARYpoint) y{i}(closestBOUNDARYpoint)];

middlePoint = repmat([1 C1 C2],size(ROI_verts,1),1);

ROI_verts_end_coords = find_point_on_line(BoundaryROIsend_coords(closest_Boundary_vert_ind,:),middlePoint,closest_Boundary_vert_dist'./max(closest_Boundary_vert_dist));

end