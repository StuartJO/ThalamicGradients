function Plot_PC123(decomp,cort_parc,savefldr)

   make_figs = 1;
if nargin < 3
   make_figs = 0; 
end

load('./data/ancillary/fsaverage_surface_data.mat')

load('./data/ancillary/MNI_Seed_voxelData.mat','seed_vox_coords','seed_mni_coords')

cmap = turbo(256);

for grad = 1:3

    pc_thal = decomp.pcs_thal(:,grad);

    PlotThalGradientSlices(pc_thal,decomp.used_seed_voxel_coords,turbo(256),['Thalamic seed PC',num2str(grad),' score'],2.1);
    
    if make_figs
    print(['./figure_outputs/',savefldr,'/Thalamus_PC',num2str(grad),'_score.png'],'-dpng','-r300')
    end

    %exportgraphics(gcf,['./figure_outputs/Thalamus_PC',num2str(grad),'_score.pdf'],'Resolution',600)

%     figure('Position',[162   233   713   592])
% 
%     SPATIAL_DIR = {'Medial-Lateral','Anterior-posterior','Dorsal-ventral'};
%     
%     for XYZ = 1:3
%         s = scatter(decomp.used_seed_mni_coords(:,XYZ),pc_thal,'filled','MarkerFaceAlpha',.25,'MarkerFaceColor',cmap(50,:));
%         [RHO,pval] = corr(decomp.used_seed_mni_coords(:,XYZ),pc_thal,'Type','Pearson');
%         xlabel({SPATIAL_DIR{XYZ},'(MNI coordinate)'})
%         ylabel(['Thalamic seed PC',num2str(grad),' score'])
%         set(gca,'FontSize',24)
%         xlimits = xlim;
%         ylimits = ylim;
%         text_x_coord = find_point_on_line(xlimits(1),xlimits(2),.05);
%         text_y_coord = find_point_on_line(ylimits(1),ylimits(2),1);
% 
%         if pval > .001
%             p_val_format = ['{\itp} = ',num2str(round(pval,3))];
%         else
%             p_val_format = '{\itp} < .001';    
%         end
%         
%         text(text_x_coord,text_y_coord,{['{\itr} = ',num2str(round(RHO,3))],p_val_format},'FontSize',24)
%         ylim([ylimits(1) find_point_on_line(ylimits(1),ylimits(2),1.1)]);
%         pval_string = num2str(pval,3);
%         
%         pval_n = pval_string(find(pval_string=='e')+1:end);
%         if make_figs
%         print(['./figure_outputs/',savefldr,'/Thalamus_PC',num2str(grad),'_',SPATIAL_DIR{XYZ},'_corr.png'],'-dpng','-r300')
%         end
%     end
%     
    
    surface.vertices = lh_inflated_verts;
    surface.faces = lh_faces;
    cmap = turbo(256);

    pc_cort = decomp.pcs_cort(:,grad);

    figure('Position',[461   462   560   325])
    ax_sub1 = axes('Position',[0.005 .33 .49 .66]);
    p = plotSurfaceROIBoundary(surface,cort_parc,pc_cort,'midpoint',cmap,2);
    camlight(80,-10);
    camlight(-80,-10);
    view([-90 0])

    axis off
    axis image

    ax_sub2 = axes('Position',[.505 .33 .489 .66]);
    plotSurfaceROIBoundary(surface,cort_parc,pc_cort,'midpoint',cmap,2);
    camlight(80,-10);
    camlight(-80,-10);
    view([90 0])
    axis off
    axis image

    caxis([min(pc_cort) max(pc_cort)])
    colormap(cmap)
    c = colorbar('Location','southoutside');
    set(c, 'Position',[.1 .23 .8 .05],'FontSize',20);
    c.Label.String = ['Cortical region PC',num2str(grad),' loading'];

    if make_figs
    print(['./figure_outputs/',savefldr,'/Cortical_PC',num2str(grad),'_coeff.png'],'-dpng','-r300')
    end

end

if make_figs
close all
end