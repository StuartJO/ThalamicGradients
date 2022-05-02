surface.vertices = lh_inflated_verts;
surface.faces = lh_faces;
vertices = lh_inflated_verts;
faces = lh_faces;
vertex_id = Scha7_parcs.lh_scha200;
% i = 74; 

% [vertices,faces] = icosphere(1);
% surface.faces = faces;
% surface.vertices = vertices;
% vertex_id = ones(size(vertices,1),1);
% vertex_id(vertices(:,1)>0) = 2;
% vertex_id(vertices(:,2)>0) = 3;


data = 1:max(vertex_id);

boundary_type = {'faces','midpoint','centroid','edges','MATLABedge','none'};
boundary_type_ind = 5;
linewidth = 2;
colorUnknownGrey = 1;

            cmap = turbo(100);

            
            figure
    [p,b,~,CHKBOUND,~,orig_data_limits] = plotSurfaceROIBoundary(surface,vertex_id,data,boundary_type{boundary_type_ind},cmap,colorUnknownGrey,linewidth);

    camlight(80,-10);
    camlight(-80,-10);

    %view([91 -10])

    axis off
    axis tight
    axis equal
%     
%         p.EdgeColor = 'k';
%     p.EdgeAlpha = .5;
    
% ylim([11.0819   17.3091])
% zlim([32.3805   34.7733])
% camva(1)

%74
p.FaceColor='flat';
hold on
%1:200

G = graph(adj);
g = plot(G);
g.XData = vertices(vert_new2old(:,2),1);g.YData = vertices(vert_new2old(:,2),2);g.ZData = vertices(vert_new2old(:,2),3);
%g.NodeCData = comp;
g.LineWidth=2;
g.EdgeColor='k';
g.EdgeAlpha=1;
g.MarkerSize=1;
%g.NodeLabel = 1:length(adj);