function PlotMouseCortFeatures(PC)

if nargin < 1
    PC = 1;
end
    
mkdir ./figure_outputs/MouseCort

mouse_CortFeatures = GetAMBAcortexdata();

load('./data/processed/mouse_decomp.mat','mouse_pcs_cort')

annotlabels = {'a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z'};
figure
for i = 1:9

   feature = mouse_CortFeatures.data(:,i); 

   %[c,p] = corr(mouse_pcs_cort(:,PC),feature,'Type','Pearson','rows','pairwise');
   
   s = scatterfit(mouse_pcs_cort(:,PC),feature,40,mouse_pcs_cort(:,PC),[],2,1);
   s.MarkerEdgeColor = [0 0 0];
   %title(['{\itr} = ',num2str(c,3),', {\itp} = ',num2str(p,3)])
   xlabel(['Cortical mouse PC',num2str(PC),' loading'])
   ylabel(mouse_CortFeatures.name{i})
   colormap(turbo(256))
   set(gca,'FontSize',16)
   a = annotation('textbox',[0 .905 .05 .13],'String',annotlabels{i},'FontSize',32,'EdgeColor','none');
   
   print(['./figure_outputs/MouseCort/mPC',num2str(PC),annotlabels{i},'.png'],'-dpng','-r300')
   close all
end