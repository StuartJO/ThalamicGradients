function Plot_PC123(decomp,cort_parc,savefldr)

   make_figs = 1;
if nargin < 3
   make_figs = 0; 
end

load('./data/ancillary/fsaverage_surface_data.mat')

cmap = turbo(256);

for grad = 1:3

    pc_thal = decomp.pcs_thal(:,grad);

    PlotThalGradientSlices(pc_thal,decomp.used_seed_voxel_coords,cmap,['Thalamic seed PC',num2str(grad),' score'],2.1);
    
    if make_figs
        print(['./figure_outputs/',savefldr,'/Thalamus_PC',num2str(grad),'_score.png'],'-dpng','-r300')
        print(['./figure_outputs/',savefldr,'/Thalamus_PC',num2str(grad),'_score.svg'],'-dsvg','-r300')
    end
    
    surface.vertices = lh_inflated_verts;
    surface.faces = lh_faces;
    pc_cort = decomp.pcs_cort(:,grad);

    figure('Position',[461   462   560*2   325*2])
    %figure('Position',[461   462   1175   638])
    ax_sub1 = axes('Position',[0.005 .33 .49 .66]);
    [p1,b1] = plotSurfaceROIBoundary(surface,cort_parc,pc_cort,'midpoint',cmap,2);
    camlight(80,-10);
    camlight(-80,-10);
    view([-90 0])

    axis off
    axis image

    ax_sub2 = axes('Position',[.505 .33 .489 .66]);
    [p2,b2] = plotSurfaceROIBoundary(surface,cort_parc,pc_cort,'midpoint',cmap,2);
    camlight(80,-10);
    camlight(-80,-10);
    view([90 0])
    axis off
    axis image

    caxis([min(pc_cort) max(pc_cort)])
    colormap(cmap)
    
    print(['./figure_outputs/',savefldr,'/Cortical_PC',num2str(grad),'_coeff_nocmap.png'],'-dpng','-r300')
    
    c = colorbar('Location','southoutside');
    set(c, 'Position',[.1 .23 .8 .05],'FontSize',38);
    c.Label.String = ['Cortical region PC',num2str(grad),' loading'];

    if make_figs
        print(['./figure_outputs/',savefldr,'/Cortical_PC',num2str(grad),'_coeff.png'],'-dpng','-r300')
        delete(p1)
        delete(p2)
        delete(b1.boundary)
        delete(b2.boundary)
        print(['./figure_outputs/',savefldr,'/Cortical_PC',num2str(grad),'_coeff.svg'],'-dsvg','-r300')
    end

end

if make_figs
    close all
end