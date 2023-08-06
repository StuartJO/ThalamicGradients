
load('./data/processed/NeuroMapCorrs.mat')
load('./data/processed/mouse_decomp.mat')
load('./data/processed/main_decomp.mat')

%%

for i = 1:3
    neuromap_corrs = load(['./data/processed/NeuroMapCorrs_PC',num2str(i),'.mat']);
    SigMapDescrips = PlotSigNeuromaps(neuromap_corrs,['./figure_outputs/SigNeuromaps_PC',num2str(i)],['Cortical region PC',num2str(i),' loading']);

    annotlabels = {'a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z'};
    NeuromapLgd{i} = '';
    for j = 1:length(SigMapDescrips)
    NeuromapLgd{i} = [NeuromapLgd{i},annotlabels{j},', ',SigMapDescrips{j},'. '];
    end
end



%% Plot PC1 vs IT

figure

decomp = load('./data/processed/main_decomp.mat');

IT_IND = find(strcmp(neuromap_corrs.name,'Intrinsic timescale'));
IT = neuromap_corrs.neuromap_parc(:,IT_IND);

pc1_cort = zscore(decomp.pcs_cort(:,1));

if neuromap_corrs.p_perm(IT_IND) == 0
     pval_format = '{\itp_{spin}} < .0001 ';   
    else
    pval_format = ['{\itp_{spin}} = ',num2str(round(neuromap_corrs.p_perm(IT_IND),4))];
end

s = scatterfit(pc1_cort,IT,40,pc1_cort,pval_format,2,1);

xlabel('Cortical region PC1 loading')
ylabel('Intrinsic timescale')
xlimits = xlim;
ylimits = ylim;
set(gca,'FontSize',18)

s.MarkerEdgeColor = [0 0 0];
colormap(turbo(256))

% text_x_coord = find_point_on_line(xlimits(1),xlimits(2),.6);
% text_y_coord = find_point_on_line(ylimits(1),ylimits(2),.9);
% text(text_x_coord,text_y_coord,{['      {\itr} = ',num2str(round(neuromap_corrs.corr(IT_IND),3))],pval_format},'FontSize',18)

DisplayCorr(pc1_cort,IT)


exportgraphics(gcf,'./figure_outputs/Fig3d.png','Resolution',300)

%%
%Plot_PC123()

decomp = load('./data/processed/main_decomp.mat');
load('./data/ancillary/fsaverage_surface_data.mat')
Plot_PC123(decomp,lh_rand500)

PlotPCs_vs_MNI(decomp,1,{'a','b','c'},'./figure_outputs/MNIcorr')
PlotPCs_vs_MNI(decomp,2,{'d','e','f'},'./figure_outputs/MNIcorr')
PlotPCs_vs_MNI(decomp,3,{'g','h','i'},'./figure_outputs/MNIcorr')

decomp = load('./data/processed/decomp_Scha400.mat');
mkdir ./figure_outputs/Scha400_PCs
Plot_PC123(decomp,Scha17_parcs.lh_scha400,'Scha400_PCs');

decomp = load('./data/processed/decomp_AllGeneSeed.mat');
mkdir ./figure_outputs/AllGeneSeed_PCs
Plot_PC123(decomp,lh_rand500,'AllGeneSeed_PCs');

%%
CompareDifferentDecomps()
%%

[MouseOhParc,MouseBrain] = GetMouseVolData;
ThalRegions = 88:122;
MouseThalOnly = MouseOhParc;
MouseThalOnly(~ismember(MouseOhParc,ThalRegions)) = NaN;
MouseThalOnly(MouseOhParc==0) = NaN;
MouseThalOnly(1:228,:,:) = NaN;

for i = 1:3
zscore_mouse_pcs_thals = zscore(mouse_score(:,i));
MouseThalPlotted = changem(MouseThalOnly,zscore_mouse_pcs_thals,ThalRegions);
PlotMouseThalGradient(MouseThalPlotted,MouseBrain,turbo(256),['Thalamic mouse PC',num2str(i),' score']);
exportgraphics(gcf,['./figure_outputs/mPC',num2str(i),'_thal.png'],'Resolution',300)
%exportgraphics(gcf,'./figure_outputs/Fig2b.png','Resolution',300)
end



%%

load('./data/ancillary/Mouseflatmap.mat', 'MouseCorticalDivisions')

for i = 1:3
figure
zscore_mouse_pcs_cort = zscore(mouse_pcs_cort(:,i));
plotData2MouseFlatmap(zscore_mouse_pcs_cort,turbo(256))

hold on

plotData2MouseFlatmap(MouseCorticalDivisions,lines(6),[1000 180 .65])

axis off
axis image

caxis([min(zscore_mouse_pcs_cort) max(zscore_mouse_pcs_cort)])
colormap(turbo(256))

c = colorbar('Location','southoutside');
c.Position = [0.4679    0.1690    0.4    0.0500];
c.Label.String = ['Cortical mouse PC',num2str(i),' loading'];

% c = colorbar('Location','southoutside');
% c.Label.String = ['Mouse cortical PC',num2str(i),' loading'];
set(gca,'FontSize',20)

exportgraphics(gcf,['./figure_outputs/MouseCortPC',num2str(i),'.png'],'Resolution',300)

end


figure
load('./data/ancillary/Mouseflatmap.mat', 'MouseCorticalDivisions')
plotData2MouseFlatmap(MouseCorticalDivisions,lines(6))
axis off
axis image
exportgraphics(gcf,'./figure_outputs/FlatMapRegions.png','Resolution',300)

%%%

TractGeneNorm = load('./data/processed/TractGeneNorm.mat');

GeneNames_human = TractGeneNorm.GeneNames;

SuppTable1_1 = table(TractGeneNorm.GeneNames,pcs_gene(:,1),pcs_gene(:,3),pcs_gene(:,2));

PhillipsGeneTable = readtable('PhillipsMouseThalGenes.xlsx');
GeneNames_phillips = PhillipsGeneTable.GeneSymbol;
GeneNames_phillips_upper = upper(GeneNames_phillips);

Phillips_pc1_gene = zscore(PhillipsGeneTable.PC1);

[~,idx] = ismember(GeneNames_human,upper(mouse_GeneNames));
AMBA2Human_homolog = idx(idx~=0);
Human2AMBA_homolog = find(idx~=0);
x = zscore(pcs_gene(Human2AMBA_homolog,1));
y = zscore(mouse_pcs_gene(AMBA2Human_homolog,1));

SuppTable1_3 = table(GeneNames_human(AMBA2Human_homolog));

figure
s = scatterfit(x,y,36,lines(1),[],2,1);
xlabel('Human gene PC1 loading')
ylabel({'Allen Mouse Brain','gene PC1 loading'})


print(['./figure_outputs/Fig2c.png'],'-dpng','-r300')

%[c,p] = corr(x,y);
%disp(['r = ',num2str(c,4),', p = ',num2str(p,5),', Human vs AMBA']);

disp('Human vs AMBA:')
DisplayCorr(x,y)

SuppTable1_2 = table(mouse_GeneNames,mouse_pcs_gene(:,1),mouse_pcs_gene(:,2),mouse_pcs_gene(:,3));

[~,idx] = ismember(GeneNames_phillips_upper,upper(mouse_GeneNames));
Phillips2AMBA_homolog = idx(idx~=0);
AMBA2Phillips_homolog = find(idx~=0);
x = zscore(Phillips_pc1_gene(AMBA2Phillips_homolog));

y = zscore(mouse_pcs_gene(Phillips2AMBA_homolog,1));

figure
s = scatterfit(x,y,36,lines(1),[],2,1);
xlabel('Phillips mouse gene PC1 loading')
ylabel({'Allen Mouse Brain','gene PC1 loading'})
print(['./figure_outputs/FigS6b.png'],'-dpng','-r300')

SuppTable1_4 = table(mouse_GeneNames(Phillips2AMBA_homolog));

% [c,p] = corr(x,y);
% disp(['r = ',num2str(c,4),', p = ',num2str(p,5),', Phillips vs AMBA']);

disp('Phillips vs AMBA:')
DisplayCorr(x,y)

[~,idx] = ismember(GeneNames_human,GeneNames_phillips_upper);
Phillips2Human_homolog = idx(idx~=0);
Human2Phillips_homolog = find(idx~=0);
x = zscore(pcs_gene(Human2Phillips_homolog,1));
y = zscore(Phillips_pc1_gene(Phillips2Human_homolog));

figure
s = scatterfit(x,y,36,lines(1),[],2,1);
ylabel({'Phillips mouse','gene PC1 loading'})
xlabel('Human gene PC1 loading')

print(['./figure_outputs/FigS6a.png'],'-dpng','-r300')
% [c,p] = corr(x,y);
% disp(['r = ',num2str(c,4),', p = ',num2str(p,5),', Human vs Phillips']);

disp('Human vs Phillips:')
DisplayCorr(x,y)

SuppTable1_5 = table(GeneNames_human(Human2Phillips_homolog));

SuppTable1_1.Properties.VariableNames = {'GeneName','PC1','PC2','PC3'};

writetable(SuppTable1_1,'SupplementaryTable1.xlsx','WriteMode','overwrite','Sheet','HumanGenesPCs')

SuppTable1_2.Properties.VariableNames = {'GeneName','mPC1','mPC2','mPC3'};

writetable(SuppTable1_2,'SupplementaryTable1.xlsx','WriteMode','Append','Sheet','MouseGenesPCs')

writetable(SuppTable1_3,'SupplementaryTable1.xlsx','WriteMode','Append','Sheet','HomologueMouseHumanGenes','WriteVariableNames',0)

writetable(SuppTable1_4,'SupplementaryTable1.xlsx','WriteMode','Append','Sheet','HomologueMousePhillipsGenes','WriteVariableNames',0)

writetable(SuppTable1_5,'SupplementaryTable1.xlsx','WriteMode','Append','Sheet','HomologueHumanPhillipsGenes','WriteVariableNames',0)


%%

for PC = 1:3
    if PC == 1
        NegativeName = 'Medial';
        PositiveName = 'Lateral';
    elseif PC == 2
        NegativeName = 'Dorsal';
        PositiveName = 'Ventral';
    elseif PC == 3
        NegativeName = 'Anterior';
        PositiveName = 'Posterior';
    end

PosTable = readtable(['./data/processed/PC',num2str(PC),'_HumanMostPositiveSpinTested.csv']);
NegTable = readtable(['./data/processed/PC',num2str(PC),'_HumanMostNegativeSpinTested.csv']);

PosTable.Properties.VariableNames = {'GeneName',['PC',num2str(PC)]};
NegTable.Properties.VariableNames = {'GeneName',['PC',num2str(PC)]};

writetable(NegTable,'SupplementaryTable2.xlsx','WriteMode','Append','Sheet',['PC',num2str(PC),'_',NegativeName,'Genes'])
writetable(PosTable,'SupplementaryTable2.xlsx','WriteMode','Append','Sheet',['PC',num2str(PC),'_',PositiveName,'Genes'])

end

%% Make 

figure

[mouse_CortFeatures,mouse_ThalHierarchy] = GetAMBAcortexdata();

s = scatterfit(mouse_pcs_cort(:,1),mouse_CortFeatures.data(:,6),40,mouse_pcs_cort(:,1),[],2,1);

disp('mPC1 hier corr:')
DisplayCorr(mouse_pcs_cort(:,1),mouse_CortFeatures.data(:,6))

xlabel('Cortical mouse PC1 loading')
ylabel('Hierarchical level')
xlimits = xlim;
ylimits = ylim;
set(gca,'FontSize',18)

s.MarkerEdgeColor = [0 0 0];
colormap(turbo(256))

exportgraphics(gcf,'./figure_outputs/Fig3e.png','Resolution',300)

close all
% 
% figure('Position',[597 266 1109 420])
% plot_label = {'A','B'};
% for i = 2:3
% subplot(1,2,i-1)
% s = scatterfit(mouse_pcs_cort(:,i),mouse_CortFeatures.data(:,6),40,mouse_pcs_cort(:,i),[],1);
% 
% xlabel(['mPC',num2str(i),' cortical loading'])
% ylabel('Hierarchical level')
% xlimits = xlim;
% ylimits = ylim;
% set(gca,'FontSize',18)
% 
% s.MarkerEdgeColor = [0 0 0];
% colormap(turbo(256))
% addPlotLabel(plot_label{i-1},gca,24,[0 0])
% end
% 
% exportgraphics(gcf,'./figure_outputs/mPCX_hier.png','Resolution',300)

close all
%%

figure
s = scatterfit(mouse_pcs_thal(:,1),mouse_ThalHierarchy,40,mouse_pcs_thal(:,1),[],2,1);
s.MarkerEdgeColor = [0 0 0];
xlabel('Thalamic mouse PC1 score')
ylabel('Thalamic hierarchy')
colormap(turbo(256))
set(gca,'FontSize',18)
%exportgraphics(gcf,'./figure_outputs/Fig2d.png','Resolution',300)
print(['./figure_outputs/Fig2d.png'],'-dpng','-r300')

disp('mPC1 thal hier corr:')
DisplayCorr(mouse_pcs_thal(:,1),mouse_ThalHierarchy)

figure('Position',[597 266 1109 420])
plot_label = {'a','b'};
for i = 2:3
subplot(1,2,i-1)
s = scatterfit(mouse_pcs_thal(:,i),mouse_ThalHierarchy,40,mouse_pcs_thal(:,i),[],2,1);

xlabel(['Thalamic mouse PC',num2str(i),' score'])
ylabel('Hierarchical level')
xlimits = xlim;
ylimits = ylim;
set(gca,'FontSize',18)

s.MarkerEdgeColor = [0 0 0];
colormap(turbo(256))
addPlotLabel(plot_label{i-1},gca,24,[0 0])
end

exportgraphics(gcf,'./figure_outputs/mPCX_hier_thal.png','Resolution',300)


%%
%plot_cell_enrichment(1)

for i = 1:3
disp(['PC',num2str(i),' cell enrichment'])
plot_cell_enrichment_PCs(i,1)
end

for i = 1:3
copyfile(['CellEnrichmentCombined_PC',num2str(i),'.xlsx'],['SupplementaryTable',num2str(i+2),'.xlsx'],'f')
end

%%
for i = 1:3
plot_gene_trajectories(i,['./figure_outputs/GeneTrajPC',num2str(i),'.png'])
end

for i = 1:3
SigDevelopmentalTrajectory = readtable(['./data/processed/SigDevelopmentalTrajectoryGenes_PC',num2str(i),'.csv']);
writetable(SigDevelopmentalTrajectory,'SupplementaryTable6.xlsx','WriteMode','Append','Sheet',['PC',num2str(i),'_','Genes'])
end

close all
%%

load('MouseThalROICoords.mat','MouseThalROICoords')
load('./data/preprocessed/AllenGeneDataset_19419.mat','structInfo')

MouseThalNuc_xcoord = MouseThalROICoords(:,1);
mouse_pc1_thal = mouse_pcs_thal(:,1);
figure('Position',[168   187   823   684])
scatter(MouseThalNuc_xcoord,mouse_pc1_thal,80,mouse_pc1_thal,'filled')
hold on
for i = 1:length(MouseThalNuc_xcoord)
    text(MouseThalNuc_xcoord((i)),mouse_pcs_thal(i,1)+.1,structInfo.acronym{(i)+87},'FontSize',12,'HorizontalAlignment','center');
end
set(gca,'FontSize',20)
xlabel('CCFv3 {\itx}-coordinate')
ylabel('Thalamic mouse PC1 score')
colormap(turbo(256))

[c,p] = corr(MouseThalNuc_xcoord,mouse_pc1_thal);
disp(['r = ',num2str(c,4),', p = ',num2str(p,5),', mouse nuclei x-position vs mPC1']);
%addPlotLabel('a')
print(['./figure_outputs/S5.svg'],'-dsvg')
print(['./figure_outputs/S5.png'],'-dpng')

PlotMousePCvsCoord(1,{'a','b','c'},'./figure_outputs/MouseCoordCorr')
PlotMousePCvsCoord(2,{'d','e','f'},'./figure_outputs/MouseCoordCorr')
PlotMousePCvsCoord(3,{'g','h','i'},'./figure_outputs/MouseCoordCorr')

for i = 1:3
    PlotMouseCortFeatures(i)
end

%%

WebGestalt_FOLDERS = {'PC1Neg','PC1Pos','PC2Neg','PC2Pos','PC3Neg','PC3Pos'};
SheetName = {'PC1_Medial','PC1_Lateral','PC2_Dorsal','PC2_Ventral','PC3_Anterior','PC3_Posterior'};
for i = 1:2
[DiseaseEnrich,DiseaseGeneTable,DiseaseEnrichResults] = ReadWebGestalt(['./data/Web/',WebGestalt_FOLDERS{i}]);
writetable(DiseaseEnrich,'SupplementaryTable7.xlsx','WriteMode','overwrite','Sheet',[SheetName{i},'DiseaseEnrich'])
writetable(DiseaseGeneTable,'SupplementaryTable7.xlsx','WriteMode','Append','Sheet',[SheetName{i},'DiseaseGenes'])
end

for i = 3:4
[DiseaseEnrich,DiseaseGeneTable,DiseaseEnrichResults] = ReadWebGestalt(['./data/Web/',WebGestalt_FOLDERS{i}]);
writetable(DiseaseEnrich,'SupplementaryTable8.xlsx','WriteMode','overwrite','Sheet',[SheetName{i},'DiseaseEnrich'])
writetable(DiseaseGeneTable,'SupplementaryTable8.xlsx','WriteMode','Append','Sheet',[SheetName{i},'DiseaseGenes'])
end

for i = 5:6
[DiseaseEnrich,DiseaseGeneTable,DiseaseEnrichResults] = ReadWebGestalt(['./data/Web/',WebGestalt_FOLDERS{i}]);
writetable(DiseaseEnrich,'SupplementaryTable9.xlsx','WriteMode','overwrite','Sheet',[SheetName{i},'DiseaseEnrich'])
writetable(DiseaseGeneTable,'SupplementaryTable9.xlsx','WriteMode','Append','Sheet',[SheetName{i},'DiseaseGenes'])
end

