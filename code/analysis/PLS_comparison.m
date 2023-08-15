% Compare the PCA result to one obtained from PLS

load('main_decomp.mat')

% Lets predict the connectivity values from the gene expression values

[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats] = plsregress(input(:,251:end),input(:,1:250));

figure('Position',[195         139        1555         459])

subplot(1,3,1)
scatter(score(:,1),XS(:,1),36,'filled')
xlabel('PC1 thalamic seed score')
ylabel({'PLS 1st component','thalamic seed score'})
set(gca,'Fontsize',16)
addPlotLabel('A',gca,24,[0 0])
subplot(1,3,2)
scatter(pcs_gene(:,1),XL(:,1),36,'filled')
xlabel('PC1 gene loading')
ylabel({'PLS 1st component','gene loading'})
set(gca,'Fontsize',16)
addPlotLabel('B',gca,24,[0 0])
subplot(1,3,3)
scatter(pcs_cort(:,1),YL(:,1),36,'filled')
xlabel('PC1 cortical loading')
ylabel({'PLS 1st component','cortical loading'})
set(gca,'Fontsize',16)
addPlotLabel('C',gca,24,[0 0])

%exportgraphics(gcf,'./figure_outputs/PLS_comparison.png')