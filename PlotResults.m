
load('./data/processed/NeuroMapCorrs.mat')
load('./data/processed/mouse_decomp.mat')
load('./data/processed/main_decomp.mat')

%%

PlotSigNeuromaps()

%% Plot PC1 vs IT

figure

IT_IND = find(strcmp(neuromap_corrs.name,'Intrinsic timescale'));
IT = neuromap_corrs.neuromap_parc(:,IT_IND);

if neuromap_corrs.p_perm(IT_IND) == 0
     pval_format = '{\itp_{spin}} < .0001 ';   
    else
    pval_format = ['{\itp_{spin}} = ',num2str(round(neuromap_corrs.p_perm(IT_IND),4))];
end

s = scatterfit(pc1_cort,IT,40,pc1_cort,pval_format,0);

xlabel('Cortical region PC1 loading')
ylabel('Intrinsic timescale')
xlimits = xlim;
ylimits = ylim;
set(gca,'FontSize',18)

s.MarkerEdgeColor = [0 0 0];
colormap(turbo(256))

text_x_coord = find_point_on_line(xlimits(1),xlimits(2),.6);
text_y_coord = find_point_on_line(ylimits(1),ylimits(2),.9);
text(text_x_coord,text_y_coord,{['      {\itr} = ',num2str(round(neuromap_corrs.corr(IT_IND),3))],pval_format},'FontSize',18)

exportgraphics(gcf,'./figure_outputs/Fig3d.png','Resolution',300)

%%
Plot_PC123()

%%
CompareDifferentDecomps()
%%

[MouseOhParc,MouseBrain] = GetMouseVolData;
ThalRegions = 88:122;
MouseThalOnly = MouseOhParc;
MouseThalOnly(~ismember(MouseOhParc,ThalRegions)) = NaN;
MouseThalOnly(MouseOhParc==0) = NaN;
MouseThalOnly(1:228,:,:) = NaN;
MouseThalPlotted = changem(MouseThalOnly,mouse_pc1_thal,ThalRegions);
PlotMouseThalGradient(MouseThalPlotted,MouseBrain,turbo(256),'mPC1 score')

exportgraphics(gcf,'./figure_outputs/Fig2b.png','Resolution',300)

%%

figure('Position',[862   208   743   630]);

plotData2MouseFlatmap(mouse_pc1_cort,turbo(256))

caxis([min(mouse_pc1_cort) max(mouse_pc1_cort)])
colormap(turbo(256))
c = colorbar('Location','southoutside');
c.Label.String = ['Mouse cortical PC1 loading'];
set(gca,'FontSize',20)

exportgraphics(gcf,'./figure_outputs/Fig3c.png','Resolution',300)

figure
load('./data/ancillary/Mouseflatmap.mat', 'MouseCorticalDivisions')
plotData2MouseFlatmap(MouseCorticalDivisions,lines(6))
exportgraphics(gcf,'./figure_outputs/FlatMapRegions.png','Resolution',300)

%%%

PhillipsGeneTable = readtable('PhillipsMouseThalGenes.xlsx');
GeneNames_phillips = PhillipsGeneTable.GeneSymbol;
GeneNames_phillips_upper = upper(GeneNames_phillips);

Phillips_pc1_gene = zscore(PhillipsGeneTable.PC1);

[~,idx] = ismember(GeneNames_human,upper(mouse_GeneNames));
AMBA2Human_homolog = idx(idx~=0);
Human2AMBA_homolog = find(idx~=0);
x = pc1_gene(Human2AMBA_homolog,1);
y = mouse_pc1_gene(AMBA2Human_homolog,1);

figure
s = scatterfit(x,y);
xlabel('Human gene PC1 loading')
ylabel({'Allen Mouse Brain','gene PC1 loading'})
print(['./figure_outputs/Fig2c.png'],'-dpng','-r300')

[c,p] = corr(x,y);
disp(['r = ',num2str(c,4),', p = ',num2str(p,5),', Human vs AMBA']);

[~,idx] = ismember(GeneNames_phillips_upper,upper(mouse_GeneNames));
Phillips2AMBA_homolog = idx(idx~=0);
AMBA2Phillips_homolog = find(idx~=0);
x = Phillips_pc1_gene(AMBA2Phillips_homolog);

y = mouse_pc1_gene(Phillips2AMBA_homolog);

figure
s = scatterfit(x,y);
xlabel('Phillips mouse gene PC1 loading')
ylabel({'Allen Mouse Brain','gene PC1 loading'})
print(['./figure_outputs/FigS6b.png'],'-dpng','-r300')

[c,p] = corr(x,y);
disp(['r = ',num2str(c,4),', p = ',num2str(p,5),', Phillips vs AMBA']);

[~,idx] = ismember(GeneNames_human,GeneNames_phillips_upper);
Phillips2Human_homolog = idx(idx~=0);
Human2Phillips_homolog = find(idx~=0);
x = pc1_gene(Human2Phillips_homolog,1);
y = Phillips_pc1_gene(Phillips2Human_homolog);

figure
s = scatterfit(x,y);
ylabel({'Phillips mouse','gene PC1 loading'})
xlabel('Human gene PC1 loading')

print(['./figure_outputs/FigS6a.png'],'-dpng','-r300')
[c,p] = corr(x,y);
disp(['r = ',num2str(c,4),', p = ',num2str(p,5),', Human vs Phillips']);


%% Make 

figure

[mouse_CortFeatures,mouse_ThalHierarchy] = GetAMBAcortexdata();

s = scatterfit(mouse_pc1_cort,mouse_CortFeatures.data(:,6),40,mouse_pc1_cort,[],1);

xlabel('Mouse cortical PC1 loading')
ylabel('Hierarchical level')
xlimits = xlim;
ylimits = ylim;
set(gca,'FontSize',18)

s.MarkerEdgeColor = [0 0 0];
colormap(turbo(256))

exportgraphics(gcf,'./figure_outputs/Fig3e.png','Resolution',300)

close all

%%

s = scatterfit(mouse_pc1_thal,mouse_ThalHierarchy,40,mouse_pc1_thal,[],1);
s.MarkerEdgeColor = [0 0 0];
xlabel('mPC1 score')
ylabel('Thalamic hierarchy')
colormap(turbo(256))
set(gca,'FontSize',18)
exportgraphics(gcf,'./figure_outputs/Fig2d.png','Resolution',300)

%%
plot_cell_enrichment(1)

%%
plot_gene_trajectories()

%%

load('MouseThalROICoord.mat','MouseThalROICoord')
load('./data/preprocessed/AllenGeneDataset_19419.mat','structInfo')

MouseThalNuc_xcoord = MouseThalROICoord(:,1);

figure('Position',[168   187   823   684])
scatter(MouseThalNuc_xcoord,mouse_pc1_thal,80,mouse_pc1_thal,'filled')
hold on
for i = 1:length(MouseThalNuc_xcoord)
    text(MouseThalNuc_xcoord((i)),mouse_pc1_thal((i))+.1,structInfo.acronym{(i)+87},'FontSize',12,'HorizontalAlignment','center');
end
set(gca,'FontSize',20)
xlabel('CCFv3 {\itx}-coordinate')
ylabel('mPC1 score')
colormap(turbo(256))

[c,p] = corr(MouseThalNuc_xcoord,mouse_pc1_thal);
disp(['r = ',num2str(c,4),', p = ',num2str(p,5),', mouse nuclei x-position vs mPC1']);
%addPlotLabel('a')
print(['./figure_outputs/S5.svg'],'-dsvg')
print(['./figure_outputs/S5.png'],'-dpng')