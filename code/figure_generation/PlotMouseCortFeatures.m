function PlotMouseCortFeatures(PC)

if nargin < 1
    PC = 1;
end
    
mkdir ./figure_outputs/MouseCort

mouse_CortFeatures = GetAMBAcortexdata();

load('./data/processed/mouse_decomp.mat','mouse_pcs_cort')

figure

for i = 1:9

   feature = mouse_CortFeatures.data(:,i); 

   annotlabel = numberToLetter(i);
     
   s = scatterfit(mouse_pcs_cort(:,PC),feature,40,mouse_pcs_cort(:,PC),[],2,1);
   s.MarkerEdgeColor = [0 0 0];
   xlabel(['Cortical mouse PC',num2str(PC),' loading'])
   ylabel(mouse_CortFeatures.name{i})
   colormap(turbo(256))
   set(gca,'FontSize',16)
   a = annotation('textbox',[0 .905 .05 .13],'String',annotlabel,'FontSize',32,'EdgeColor','none');
      
   [RHO,P,CIL,CIU] = corrcoef(mouse_pcs_cort(:,PC),feature,'Rows','complete');
   
   DF = sum(~isnan(sum(data,2)))-2;
    
   disp(['Panel ',annotlabel,', Pearson''s r(',num2str(DF),') = ',num2str(RHO(1,2)),', p = ',num2str(P(1,2)),', CI = [',num2str(CIL(1,2)),', ',num2str(CIU(1,2)),'], two-tailed'])
      
   print(['./figure_outputs/MouseCort/mPC',num2str(PC),annotlabel,'.png'],'-dpng','-r300')
   print(['./figure_outputs/MouseCort/mPC',num2str(PC),annotlabel,'.svg'],'-dsvg')
   close all
end