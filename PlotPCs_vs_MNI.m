function PlotPCs_vs_MNI(decomp,grad,annotlabels,SavePreFix)

SPATIAL_DIR = {'Medial-lateral','Anterior-posterior','Dorsal-ventral'};
SPATIAL_AXIS = {'x','y','z'};

pc_thal = decomp.pcs_thal(:,grad);

for XYZ = 1:3
    figure('Position',[162   233   713   592])
    
    s = scatterfit(decomp.used_seed_mni_coords(:,XYZ)*-1,pc_thal,36,[0.2542    0.5895    0.9990],[],0);
    s.MarkerFaceAlpha = .25;
    %s = scatter(decomp.used_seed_mni_coords(:,XYZ)*-1,pc_thal,'filled','MarkerFaceAlpha',.25,'MarkerFaceColor',[0.2542    0.5895    0.9990]);
    [R,pval] = corr(decomp.used_seed_mni_coords(:,XYZ)*-1,pc_thal,'Type','Pearson');
    
    [RHO,RHOpval] = corr(decomp.used_seed_mni_coords(:,XYZ)*-1,pc_thal,'Type','Spearman');
    
    [~,~,CIL,CIU] = corrcoef(decomp.used_seed_mni_coords(:,XYZ)*-1,pc_thal);
    
    FigCaption{XYZ} = [annotlabels{XYZ},', PC',num2str(grad),' score with ',lower(SPATIAL_DIR{XYZ}),' axis (MNI ',SPATIAL_AXIS{XYZ},'-coordinate): Pearson''s r(',num2str(length(pc_thal)-1),') = ',num2str(R),', p = ',num2str(pval),', CI = [',num2str(CIL(1,2)),', ',num2str(CIU(1,2)),'], two-tailed.'];
             
    disp(FigCaption{XYZ})
    
    SpearmanResult = ['PC',num2str(grad),' score with ',lower(SPATIAL_DIR{XYZ}),' axis (MNI ',SPATIAL_AXIS{XYZ},'-coordinate): Spearmans''s rs(',num2str(length(pc_thal)-1),') = ',num2str(RHO),', p = ',num2str(RHOpval),', two-tailed.'];
    disp(SpearmanResult)
    
    xlabel({SPATIAL_DIR{XYZ},['(MNI {\it',SPATIAL_AXIS{XYZ},'}-coordinate)']})
    ylabel(['Thalamic seed PC',num2str(grad),' score'])
    set(gca,'FontSize',24)
    xlimits = xlim;
    ylimits = ylim;

    if pval >= .001
        p_val_format = ['{\itp} = ',num2str(round(pval,3))];
    else
        %p_val_format = '{\itpval} < .001';
        p_rounded = num2str(pval,'%.2s');        
        eloc = strfind(p_rounded,'e');    
        p_val_format = ['{\itp} = ',p_rounded(1:eloc-1),'\times10^{',p_rounded(eloc+1:end),'}'];            
    end
    
    title(['{\itr} = ',num2str(round(R,3)),', ',p_val_format],'FontWeight','Normal')
    a = annotation('textbox',[0 .895 .05 .13],'String',annotlabels{XYZ},'FontSize',32,'EdgeColor','none');
    
    print([SavePreFix,'_',annotlabels{XYZ},'.png'],'-dpng','-r300')
    
end