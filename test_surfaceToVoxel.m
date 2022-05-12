coeff{3}(1:250,1)

x = repmat(1:2,250,1);
y = [(1:250)' (1:250)'];

plot(y)


verts_dist2centre = zeros(length(lh_rand500),1);
for i = 1:250
    roi_verts_ind = find(lh_rand500==i);
    roi_verts = surface.vertices(roi_verts_ind,:);
    d = squareform(pdist(roi_verts));
    [~,I] = min(mean(d));
    ROI_centre(i,:) = roi_verts(I,:);
    verts_dist2centre(roi_verts_ind) = d(I,:);
end

p = plotSurfaceROIBoundary(surface,lh_rand500,1:250,'midpoint',turbo(250),1,2);

surface2 = surface;
surface2.vertices = ras';

p2 = plotSurfaceROIBoundary(surface2,lh_rand500,1:250,'midpoint',turbo(250),1,2);


surface2.vertices = lh_mni_vox;
[p2,b,BOUNDARY] = plotSurfaceROIBoundary(surface2,lh_rand500,1:250,'midpoint',turbo(250),2);
hold on
h = imagesc(brain(:,:,90)');
colormap(gray(256))
camlight(80,-10);
camlight(-80,-10);
view([90 0])
%axis off
axis tight
axis equal
axis vis3d
p2.Vertices(:,3) = lh_mni_vox(:,3)-90;
for i = 1:251
    b.boundary(i).ZData = BOUNDARY{i}(:,3)-90;
end


surface2.vertices = lh_mni_vox;
[p2,b,BOUNDARY] = plotSurfaceROIBoundary(surface2,lh_rand500,1:250,'midpoint',turbo(250),2);

camlight(80,-10);
camlight(-80,-10);
view([90 0])
xlimits = xlim;
ylimits = ylim;
zlimits = zlim;
%axis off
axis tight
axis equal
axis vis3d
hold on
brain_slice = brain(:,:,90)';
colormap(gray(256))
slicesurf = surf(repmat(z_loc, size(brain_slice)), brain_slice);
slicesurf.EdgeAlpha=0;
colormap(gray(256))
xlim(xlimits)
ylim(ylimits)
zlim(zlimits)
axis tight