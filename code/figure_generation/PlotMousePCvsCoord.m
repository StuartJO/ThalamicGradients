function FigCaption = PlotMousePCvsCoord(grad,annotlabels,SavePreFix)

%annotlabels = {'a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z'};

SPATIAL_AXIS = {'x','y','z'};
SPATIAL_DIR = {'Medial-lateral','Anterior-posterior','Ventral-dorsal'};
cmap = turbo(256);

load('./data/processed/mouse_decomp.mat','mouse_pcs_thal')
load('MouseThalROICoords.mat','MouseThalROICoords')

FigCaption = cell(3,1);

iter = 1; 
    mpc_thal = zscore(mouse_pcs_thal(:,grad));
    for XYZ = 1:3
        %subplot(3,3,sub2ind([3 3],XYZ,grad))
        figure('Position',[162   233   713   592])
            %s = scatter(MouseThalROICoords(:,XYZ),mpc_thal,40,'filled','MarkerFaceAlpha',1,'MarkerFaceColor',cmap(50,:));
            
            s = scatterfit(MouseThalROICoords(:,XYZ),mpc_thal,80,cmap(50,:),[],0);
            
            [RHO,pval] = corr(MouseThalROICoords(:,XYZ),mpc_thal,'Type','Pearson');
            
            [~,~,CIL,CIU] = corrcoef(MouseThalROICoords(:,XYZ),mpc_thal);
    
            FigCaption{XYZ} = [annotlabels{XYZ},', mPC',num2str(grad),' score with ',lower(SPATIAL_DIR{XYZ}),' axis (CCFv3 ',SPATIAL_AXIS{XYZ},'-coordinate): Pearson''s r(',num2str(length(mpc_thal)-2),') = ',num2str(RHO),', p = ',num2str(pval),', CI = [',num2str(CIL(1,2)),', ',num2str(CIU(1,2)),'], two-tailed.'];
             
            disp(FigCaption{XYZ})

            xlabel({SPATIAL_DIR{XYZ},['(CCFv3 {\it',SPATIAL_AXIS{XYZ},'}-coordinate)']})
            ylabel(['Thalamic mouse PC',num2str(grad),' score'])
            set(gca,'FontSize',24)
            xlimits = xlim;
            ylimits = ylim;
            text_x_coord = find_point_on_line(xlimits(1),xlimits(2),.05);
            text_y_coord = find_point_on_line(ylimits(1),ylimits(2),1);

            if pval > .001
                p_val_format = ['{\itp} = ',num2str(round(pval,3))];
            else
                %p_val_format = '{\itp} < .001';  
                 p_rounded = num2str(pval,'%.2s');        
                 eloc = strfind(p_rounded,'e');    
                 p_val_format = ['{\itp} = ',p_rounded(1:eloc-1),'\times10^{',p_rounded(eloc+1:end),'}'];            
            end
            title(['{\itr} = ',num2str(round(RHO,3)),', ',p_val_format],'FontWeight','Normal')
            a = annotation('textbox',[0 .896 .05 .13],'String',annotlabels{XYZ},'FontSize',32,'EdgeColor','none');
            print([SavePreFix,'_',annotlabels{XYZ},'.png'],'-dpng','-r300')
            print([SavePreFix,'_',annotlabels{XYZ},'.svg'],'-dsvg')
    end

end