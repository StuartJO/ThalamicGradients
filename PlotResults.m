
%% This script makes all the figures


load('./data/processed/mouse_decomp.mat')
load('./data/processed/decomp_rand500.mat')
%load('./data/preprocessed/AllenGeneDataset_19419.mat','structInfo')

mkdir ./SourceDataTables
mkdir ./SupplementaryTables
mkdir ./figure_outputs
%warning('off',w.identifier)

%% Plot PC1, 2, and 3 for different decompositions (Figs: 2a, 3a, S2, S6, S15

Seeds = (1:length(pcs_thal(:,1)))';

CorticalROIs = (1:length(pcs_cort(:,1)))';

disp('Plot human PC results')

decomp_rand500 = load('./data/processed/decomp_rand500.mat');
load('./data/ancillary/fsaverage_surface_data.mat')
mkdir ./figure_outputs/Rand500_PCs
Plot_PC123(decomp_rand500,lh_rand500,'Rand500_PCs')

mkdir ./figure_outputs/MNIcorr
PlotPCs_vs_MNI(decomp_rand500,1,{'a','b','c'},'./figure_outputs/MNIcorr/');
PlotPCs_vs_MNI(decomp_rand500,2,{'d','e','f'},'./figure_outputs/MNIcorr/');
PlotPCs_vs_MNI(decomp_rand500,3,{'g','h','i'},'./figure_outputs/MNIcorr/');

TableData = [decomp_rand500.pcs_thal(:,1:3) (decomp_rand500.used_seed_mni_coords*-1)];
SourceDataTable = array2table(TableData,'VariableNames',{'PC1_score','PC2_score','PC3_score','MNI_x_coord','MNI_y_coord','MNI_z_coord'});
writetable(SourceDataTable,'./SourceDataTables/FigS1.xlsx')

decompScha400 = load('./data/processed/decomp_Scha400.mat');
mkdir ./figure_outputs/Scha400_PCs
Plot_PC123(decompScha400,Scha17_parcs.lh_scha400,'Scha400_PCs');

SourceDataTable = array2table([(1:921)' decompScha400.pcs_thal(:,1)],'VariableNames',{'Thalamic_Seed','PC1_score'});
writetable(SourceDataTable,'./SourceDataTables/FigS15b.xlsx')
SourceDataTable = array2table([(1:200)' decompScha400.pcs_cort(:,1)],'VariableNames',{'Cortical_ROI','PC1_score'});
writetable(SourceDataTable,'./SourceDataTables/FigS15a.xlsx')

decompAllGene = load('./data/processed/decomp_AllGeneSeed.mat');
mkdir ./figure_outputs/AllGeneSeed_PCs
Plot_PC123(decompAllGene,lh_rand500,'AllGeneSeed_PCs');

SourceDataTable = array2table([ [1:length(decompAllGene.pcs_thal(:,1))]' decompAllGene.pcs_thal(:,1)],'VariableNames',{'Thalamic_Seed','PC1_score'});
writetable(SourceDataTable,'./SourceDataTables/FigS6a.xlsx')
SourceDataTable = array2table([CorticalROIs decompAllGene.pcs_cort(:,1)],'VariableNames',{'Cortical_ROI','PC1_score'});
writetable(SourceDataTable,'./SourceDataTables/FigS6b.xlsx')

%% Make Mouse PC1, 2 and 3 (Figs 2b, S11a, S11b)

disp('Plot mouse thalamic gradients')

[MouseOhParc,MouseBrain] = GetMouseVolData;
ThalRegions = 88:122;
MouseThalOnly = MouseOhParc;
MouseThalOnly(~ismember(MouseOhParc,ThalRegions)) = NaN;
MouseThalOnly(MouseOhParc==0) = NaN;
MouseThalOnly(1:228,:,:) = NaN;

mkdir './figure_outputs/mPCs'

for i = 1:3
    zscore_mouse_pcs_thals = zscore(mouse_score(:,i));
    MouseThalPlotted = changem(MouseThalOnly,zscore_mouse_pcs_thals,ThalRegions);
    PlotMouseThalGradient(MouseThalPlotted,MouseBrain,turbo(256),['Thalamic mouse PC',num2str(i),' score']);
    exportgraphics(gcf,['./figure_outputs/mPCs/mPC',num2str(i),'_thal.png'],'Resolution',300)
    print(['./figure_outputs/mPCs/mPC',num2str(i),'_thal.svg'],'-dsvg','-r300')
        
end

%% Compare human, mouse, and Phillips genes (Figs 2c, S12)

disp('Plot comparison of human, mouse, and Phillips genes')

TractGeneNorm = load('./data/processed/TractGeneNorm_rand500.mat');

GeneNames_human = TractGeneNorm.GeneNames_human;

SuppTable5_1 = table(TractGeneNorm.GeneNames_human,pcs_gene(:,1),pcs_gene(:,3),pcs_gene(:,2));

PhillipsGeneTable = readtable('./data/preprocessed/PhillipsMouseThalGenes.xlsx');
GeneNames_phillips = PhillipsGeneTable.GeneSymbol;
GeneNames_phillips_upper = upper(GeneNames_phillips);

Phillips_pc1_gene = zscore(PhillipsGeneTable.PC1);

% Compare human and AMBA gene PC1 loadings

[~,idx] = ismember(GeneNames_human,upper(GeneNames_mouse));
AMBA2Human_homolog = idx(idx~=0);
%Human2AMBA_homolog = find(idx~=0);
Human2AMBA_homolog = idx~=0;
x = zscore(pcs_gene(Human2AMBA_homolog,1));
y = zscore(mouse_pcs_gene(AMBA2Human_homolog,1));

SourceDataTable = cell2table([GeneNames_human(AMBA2Human_homolog) num2cell([x y])],'VariableNames',{'Gene','Human_PC1_loading','AMBA_PC1_loading'},'RowNames',GeneNames_human(AMBA2Human_homolog));
writetable(SourceDataTable,'./SourceDataTables/Fig2c.xlsx')

SuppTable5_3 = table(GeneNames_human(AMBA2Human_homolog));

figure
s = scatterfit(x,y,36,lines(1),[],2,1);
xlabel('Human gene PC1 loading')
ylabel({'Allen Mouse Brain','gene PC1 loading'})

print('./figure_outputs/Fig2c.png','-dpng','-r300')
print('./figure_outputs/Fig2c.svg','-dsvg','-r300')

disp('Human vs AMBA:')
DisplayCorr(x,y)

SuppTable5_2 = table(GeneNames_mouse,mouse_pcs_gene(:,1),mouse_pcs_gene(:,2),mouse_pcs_gene(:,3));

% Compare Phillips and AMBA gene PC1 loadings

[~,idx] = ismember(GeneNames_phillips_upper,upper(GeneNames_mouse));
Phillips2AMBA_homolog = idx(idx~=0);
AMBA2Phillips_homolog = find(idx~=0);
x = zscore(Phillips_pc1_gene(AMBA2Phillips_homolog));

y = zscore(mouse_pcs_gene(Phillips2AMBA_homolog,1));

figure
s = scatterfit(x,y,36,lines(1),[],2,1);
xlabel('Phillips mouse gene PC1 loading')
ylabel({'Allen Mouse Brain','gene PC1 loading'})
print('./figure_outputs/FigS12b.png','-dpng','-r300')
print('./figure_outputs/FigS12b.svg','-dsvg','-r300')

SuppTable5_4 = table(GeneNames_mouse(Phillips2AMBA_homolog));

disp('Phillips vs AMBA:')
DisplayCorr(x,y)
SourceDataTable = cell2table([GeneNames_mouse(Phillips2AMBA_homolog) num2cell([x y])],'VariableNames',{'Gene','Phillips_PC1_loading','AMBA_PC1_loading'},'RowNames',GeneNames_mouse(Phillips2AMBA_homolog));
writetable(SourceDataTable,'./SourceDataTables/FigS12b.xlsx')

% Compare human and Phillips gene PC1 loadings

[~,idx] = ismember(GeneNames_human,GeneNames_phillips_upper);
Phillips2Human_homolog = idx(idx~=0);
Human2Phillips_homolog = find(idx~=0);
x = zscore(pcs_gene(Human2Phillips_homolog,1));
y = zscore(Phillips_pc1_gene(Phillips2Human_homolog));

figure
s = scatterfit(x,y,36,lines(1),[],2,1);
ylabel({'Phillips mouse','gene PC1 loading'})
xlabel('Human gene PC1 loading')

print('./figure_outputs/FigS12a.png','-dpng','-r300')
print('./figure_outputs/FigS12a.svg','-dsvg','-r300')

disp('Human vs Phillips:')
DisplayCorr(x,y)

SourceDataTable = cell2table([GeneNames_human(Human2Phillips_homolog) num2cell([x y])],'VariableNames',{'Gene','Human_PC1_loading','Phillips_PC1_loading'},'RowNames',GeneNames_human(Human2Phillips_homolog));
writetable(SourceDataTable,'./SourceDataTables/FigS12a.xlsx')

SuppTable5_5 = table(GeneNames_human(Human2Phillips_homolog));

SuppTable5_1.Properties.VariableNames = {'GeneName','PC1','PC2','PC3'};

filename = './SupplementaryTables/SupplementaryTable5.xlsx';
if exist(filename, 'file')
    delete(filename);
    disp(['File "', filename, '" deleted. Remaking.']);
else
    disp(['Making "', filename,'"']);
end

writetable(SuppTable5_1,'./SupplementaryTables/SupplementaryTable5.xlsx','WriteMode','overwrite','Sheet','HumanGenesPCs')

SuppTable5_2.Properties.VariableNames = {'GeneName','mPC1','mPC2','mPC3'};

writetable(SuppTable5_2,'./SupplementaryTables/SupplementaryTable5.xlsx','WriteMode','Append','Sheet','MouseGenesPCs')

writetable(SuppTable5_3,'./SupplementaryTables/SupplementaryTable5.xlsx','WriteMode','Append','Sheet','HomologueMouseHumanGenes','WriteVariableNames',0)

writetable(SuppTable5_4,'./SupplementaryTables/SupplementaryTable5.xlsx','WriteMode','Append','Sheet','HomologueMousePhillipsGenes','WriteVariableNames',0)

writetable(SuppTable5_5,'./SupplementaryTables/SupplementaryTable5.xlsx','WriteMode','Append','Sheet','HomologueHumanPhillipsGenes','WriteVariableNames',0)

%% Make Fig 2d, S14

disp('Make plot of mouse thalamic hierarchy vs PCs')

figure

[~,mouse_ThalHierarchy] = GetAMBAcortexdata();

s = scatterfit(mouse_pcs_thal(:,1),mouse_ThalHierarchy,40,mouse_pcs_thal(:,1),[],2,1);
s.MarkerEdgeColor = [0 0 0];
xlabel('Thalamic mouse PC1 score')
ylabel('Thalamic hierarchy')
colormap(turbo(256))
set(gca,'FontSize',18)
print('./figure_outputs/Fig2d.png','-dpng','-r300')
print('./figure_outputs/Fig2d.svg','-dsvg','-r300')

SourceDataTable = cell2table([mThalROIs num2cell([mouse_pcs_thal(:,1) mouse_ThalHierarchy])],'VariableNames',{'Thalamic_ROI','mPC1_score','Hierarchical_level'},'RowNames',mThalROIs);
writetable(SourceDataTable,'./SourceDataTables/Fig2d.xlsx')

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

disp('mPC2 thal hier corr:')
DisplayCorr(mouse_pcs_thal(:,2),mouse_ThalHierarchy)

disp('mPC3 thal hier corr:')
DisplayCorr(mouse_pcs_thal(:,3),mouse_ThalHierarchy)

exportgraphics(gcf,'./figure_outputs/FigS14.png','Resolution',300)
print('./figure_outputs/FigS14.svg','-dsvg','-r300')

SourceDataTable = cell2table([mThalROIs num2cell([mouse_pcs_thal(:,2) mouse_ThalHierarchy])],'VariableNames',{'Thalamic_ROI','mPC2_score','Hierarchical_level'},'RowNames',mThalROIs);
writetable(SourceDataTable,'./SourceDataTables/FigS14a.xlsx')

SourceDataTable = cell2table([mThalROIs num2cell([mouse_pcs_thal(:,3) mouse_ThalHierarchy])],'VariableNames',{'Thalamic_ROI','mPC3_score','Hierarchical_level'},'RowNames',mThalROIs);
writetable(SourceDataTable,'./SourceDataTables/FigS14b.xlsx')

%% Make Fig 3b

PlotTractogramPCs()

%% Make the mouse flatmap plots (Figs 3c, S11c, S11d)

disp('Plot mouse cortical gradients')

load('./data/ancillary/Mouseflatmap.mat', 'MouseCorticalDivisions')

for i = 1:3
    figure
    zscore_mouse_pcs_cort = zscore(mouse_pcs_cort(:,i));
    FlatMapPC = plotData2MouseFlatmap(zscore_mouse_pcs_cort,turbo(256));

    hold on

    FlatMapROIs = plotData2MouseFlatmap(MouseCorticalDivisions,lines(6),[1000 180 .65]);

    axis off
    axis image
    exportgraphics(gcf,['./figure_outputs/mPCs/MouseCortPC',num2str(i),'_no_cbar.png'],'Resolution',1200)
    caxis([min(zscore_mouse_pcs_cort) max(zscore_mouse_pcs_cort)])
    colormap(turbo(256))

    c = colorbar('Location','southoutside');
    c.Position = [0.4679    0.1690    0.4    0.0500];
    c.Label.String = ['Cortical mouse PC',num2str(i),' loading'];

    set(gca,'FontSize',20)

    exportgraphics(gcf,['./figure_outputs/mPCs/MouseCortPC',num2str(i),'.png'],'Resolution',1200)
    delete(FlatMapROIs.plot)
    delete(FlatMapPC.plot)
    print(['./figure_outputs/MouseCortPC',num2str(i),'_cbar.svg'],'-dsvg','-r300')
end

figure
load('./data/ancillary/Mouseflatmap.mat', 'MouseCorticalDivisions')
plotData2MouseFlatmap(MouseCorticalDivisions,lines(6));
axis off
axis image
exportgraphics(gcf,'./figure_outputs/mPCs/FlatMapRegions.png','Resolution',300)

%% Plot Figure 3d (PC1 vs IT)

disp('Getting PC1 vs IT')

figure

decomp_rand500 = load('./data/processed/decomp_rand500.mat');
neuromap_corrs = load('./data/processed/NeuroMapCorrs_PC1.mat');

IT_IND = find(strcmp(neuromap_corrs.name,'Intrinsic timescale'));
IT = neuromap_corrs.neuromap_parc(:,IT_IND);

pc1_cort = zscore(decomp_rand500.pcs_cort(:,1));

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

DisplayCorr(pc1_cort,IT)

exportgraphics(gcf,'./figure_outputs/Fig3d.png','Resolution',300)
print('./figure_outputs/Fig3d.svg','-dsvg','-r300')

SourceDataTable = array2table([pc1_cort IT],'VariableNames',{'PC1_loading','Intrinsic_timescale'});
   
writetable(SourceDataTable,'./SourceDataTables/Fig3d.xlsx')

%% Make Figure 3e 

disp('Make plot of mouse cortical hierarchy vs PC1')

figure

mouse_CortFeatures = GetAMBAcortexdata();

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
print('./figure_outputs/Fig3e.svg','-dsvg','-r300')

SourceDataTable = cell2table([mCortROIs num2cell([mouse_pcs_cort(:,1) mouse_CortFeatures.data(:,6)])],'VariableNames',{'Cortical_ROI','mPC1_loading','Hierarchical_level'},'RowNames',mCortROIs);
writetable(SourceDataTable,'./SourceDataTables/Fig3e.xlsx')

close all

%% Make cell enrichment plots (Figs 4, S22, S23)


disp('Make cell enrichment plots')

for i = 1:3
disp(['PC',num2str(i),' cell enrichment'])
plot_cell_enrichment_PCs(i,['./data/processed/CellEnrichmentCombined_PC',num2str(i),'.xlsx'])
exportgraphics(gcf,['./figure_outputs/cell_enrich_PC',num2str(i),'.png'],'Resolution',300)
print(['./figure_outputs/cell_enrich_PC',num2str(i),'.svg'],'-dsvg')

% The above output is massive, because it is rendering every single cell in
% the tsne plot. To fix this, we just remove those parts of the file

deleteSvgLines(['./figure_outputs/cell_enrich_PC',num2str(i),'.svg']);

% However we still want to show the figure. We can use the .png image as
% this uses up less memory/space. What we do is delete everything from the
% figure that isn't the tsne plot. The saved output of this figure should
% have the exact same dimensions as the .svg made above, so they can be
% easily combined in seperate software so you will still have pretty vector
% graphics for all the non-tsne parts, with an image of the tsne plot
% instead

f1 = gcf;
allAxesInFigure = findall(f1,'type','axes');
delete(allAxesInFigure(7:9))

allTextBoxInFigure = findall(f1,'type','textbox');
delete(allTextBoxInFigure)

for j = 1:6
    delete(findall(allAxesInFigure(j),'type','text'))
end

print(['./figure_outputs/cell_enrich_PC',num2str(i),'_tsne.png'],'-dpng','-r600')

end

FigNumber = {'4','S22','S23'};

tsnedata = readtable('tsne_subdata.xlsx');

SubClustID_char = tsnedata.cx;
ClusterID_SubID = cell2mat(cellfun(@str2num,(split(SubClustID_char,'-')),'UniformOutput',false));

TSNENeuronSubclustersIND = ClusterID_SubID(:,1)<=3;

tsne_neuron_only = tsnedata(TSNENeuronSubclustersIND,:);

NeuronSubtype = {'Habenula','Rora','Gad2/Ahi1'};

SuppTableNums = 7:9;

for i = 1:3
    filename = ['./SupplementaryTables/SupplementaryTable',num2str(SuppTableNums(i)),'.xlsx'];
    if exist(filename, 'file')
        delete(filename);
        disp(['File "', filename, '" deleted. Remaking.']);
    else
        disp(['Making "', filename,'"']);
    end
copyfile(['./data/processed/CellEnrichmentCombined_PC',num2str(i),'.xlsx'],['./SupplementaryTables/SupplementaryTable',num2str(SuppTableNums(i)),'.xlsx'],'f')

myTable = readtable(['CellEnrichmentCombined_PC',num2str(i),'.xlsx'],'Sheet','CellClassEnrichment');
writetable(myTable,['./SourceDataTables/Fig',FigNumber{i},'a.xlsx'])

myTable = readtable(['CellEnrichmentCombined_PC',num2str(i),'.xlsx'],'Sheet','NeuronSubtypeEnrichment');
writetable(myTable,['./SourceDataTables/Fig',FigNumber{i},'b.xlsx'])

myTable = readtable(['CellEnrichmentCombined_PC',num2str(i),'.xlsx'],'Sheet','NeuronSubclusterEnrichment');
writetable(myTable,['./SourceDataTables/Fig',FigNumber{i},'c.xlsx'])

Neuronal_subclusterResults = table2array(myTable(:,2:5));

NeuronSubtypeID = cell(length(myTable.Neuronal_subcluster),1);

tsne_neuron_enrichment = zeros(sum(TSNENeuronSubclustersIND),2);

    for j = 1:length(myTable.Neuronal_subcluster)    
        if contains(myTable.Neuronal_subcluster{j},NeuronSubtype{1})
            NeuronSubtypeID{j} = strrep(myTable.Neuronal_subcluster{j},NeuronSubtype{1},'1');
        elseif contains(myTable.Neuronal_subcluster{j},NeuronSubtype{2})
            NeuronSubtypeID{j} = strrep(myTable.Neuronal_subcluster{j},NeuronSubtype{2},'2');
        elseif contains(myTable.Neuronal_subcluster{j},NeuronSubtype{3})
            NeuronSubtypeID{j} = strrep(myTable.Neuronal_subcluster{j},NeuronSubtype{3},'3');
        end   
        
        IND = contains(tsne_neuron_only.cx,NeuronSubtypeID{j});
        tsne_neuron_enrichment(IND,1) = Neuronal_subclusterResults(j,1);
        tsne_neuron_enrichment(IND,2) = Neuronal_subclusterResults(j,3);
    end

    SourceDataTable = array2table([tsne_neuron_only.V1 tsne_neuron_only.V2 tsne_neuron_enrichment],'VariableNames',{'TSNE_x_coord','TSNE_y_coord',myTable.Properties.VariableNames{2},myTable.Properties.VariableNames{4}});
    writetable(SourceDataTable,['./SourceDataTables/Fig',FigNumber{i},'d.xlsx'])
       
end

close all

%% Plot the alternative decompositions (Figs S3, S4)

disp('Plot alternative decomposition')

CompareDifferentDecomps()

load('./data/processed/alt_decomps.mat','AltEmbedding')

Tabledata = zeros(length(Seeds),length(AltEmbedding));

ResultsTableHeader1 = cell(length(AltEmbedding),1);
ResultsTableHeader2 = cell(length(AltEmbedding),1);

for i = 1:length(AltEmbedding)
    Tabledata(:,i) = AltEmbedding{i}(:,1);
    ResultsTableHeader1{i} = ['FigS3',numberToLetter(i)]; 
    ResultsTableHeader2{i} = ['FigS4',numberToLetter(i),'_ydata']; 
end

ResultsTableHeader = [{'Thalamic_Seed'}; ResultsTableHeader1];
SourceDataTable = array2table([Seeds Tabledata],'VariableNames',ResultsTableHeader);
writetable(SourceDataTable,'./SourceDataTables/FigS3.xlsx')

ResultsTableHeader = [{'FigS4_xdata'}; ResultsTableHeader2];
SourceDataTable = array2table([Seeds Tabledata],'VariableNames',ResultsTableHeader);
writetable(SourceDataTable,'./SourceDataTables/FigS4.xlsx')


%% Make Fig. S5

MainGeneSeed = load('./data/processed/TractGeneNorm_rand500.mat');
AllGeneSeed = load('./data/processed/TractGeneNorm_AllGeneSeed.mat');

load('./data/preprocessed/MNI_Seed_voxelData.mat','seed_vox_coords','seed_mni_coords')

seedIn = MainGeneSeed.seed_ind + AllGeneSeed.seed_ind;

CoverageCmap = [1 1 1; lines(2)];

c = PlotThalGradientSlices(seedIn,seed_vox_coords,CoverageCmap,'Coverage',2.1);

caxis([-.5 2.5])

c.Ticks = 0:2;
c.TickLabels = {'Full','Gene','Used'};

print('./figure_outputs/ThalamicCoverage','-dpng','-r300')
print('./figure_outputs/FigureS5.svg','-dsvg','-r300')

SourceDataTable = array2table([seed_mni_coords AllGeneSeed.seed_ind' MainGeneSeed.seed_ind'],'VariableNames',{'Seed_MNI_x_coord','Seed_MNI_y_coord','Seed_MNI_z_coord','Seed_in_gene_coverage','Seed_in_main_analysis'});
writetable(SourceDataTable,'./SourceDataTables/FigS5.xlsx')

close all

%% Plot bootstrap result (Fig S7)

tract_decomp = load('decomp_TractOnly.mat');

decomp_rand500 = load('decomp_rand500.mat');

load('Bootstrap_results.mat');

figure('Position',[146 301 1073 590])
offset_x = -.01;
offset_y = .02;

subplot(2,3,2)
histogram(corr(decomp_rand500.score(:,1),scorePC1_boot))
xlabel('Correlation with PC1 score')
ylabel('Bootstrap count')
title('Gene+connectivity')
addPlotLabel('b',gca,24,[offset_x offset_y])
subplot(2,3,3)
histogram(corr(decomp_rand500.coeff(:,1),coeffPC1_boot))
xlabel('Correlation with PC1 loading')
ylabel('Bootstrap count')
title('Gene+connectivity')
addPlotLabel('c',gca,24,[offset_x offset_y])
subplot(2,3,5)
histogram(corr(tract_decomp.score(:,1),scorePC1_boot_tract))
xlabel('Correlation with PC1 score')
ylabel('Bootstrap count')
title('Connectivity')
addPlotLabel('e',gca,24,[offset_x offset_y])
subplot(2,3,6)
histogram(corr(tract_decomp.coeff(:,1),coeffPC1_boot_tract))
xlabel('Correlation with PC1 loading')
ylabel('Bootstrap count')
title('Connectivity')
addPlotLabel('f',gca,24,[offset_x offset_y])

subplot(2,3,1)
for i = 1:5
V(i) = Violin(expl_boot(:,i), i,'ShowMean',true,'ShowData',false,'ViolinAlpha',.25);
V(i).ViolinColor = [.5 .5 .5];
V(i).MeanPlot.Color = [0 0 0];
V(i).EdgeColor = [0 0 0];
V(i).ViolinPlot.LineWidth = .5;
V(i).BoxColor = [0 0 0];
end
hold on
plot(decomp_rand500.explained(1:5),'r','LineWidth',2)
xlabel('PC')
ylabel('Variance explained')
title('Gene+connectivity')
addPlotLabel('a',gca,24,[offset_x offset_y])
xticks(1:5)
subplot(2,3,4)
for i = 1:5
V(i) = Violin(expl_boot_tract(:,i), i,'ShowMean',true,'ShowData',false,'ViolinAlpha',.25);
V(i).ViolinColor = [.5 .5 .5];
V(i).MeanPlot.Color = [0 0 0];
V(i).EdgeColor = [0 0 0];
V(i).ViolinPlot.LineWidth = .5;
V(i).BoxColor = [0 0 0];
end
hold on
plot(tract_decomp.explained(1:5),'r','LineWidth',2)
xlabel('PC')
title('Connectivity')
ylabel('Variance explained')
addPlotLabel('d',gca,24,[offset_x offset_y])
xticks(1:5)

exportgraphics(gcf,'./figure_outputs/FigureS7.png','Resolution',300)
print('./figure_outputs/FigureS7.svg','-dsvg')

Nboots = size(expl_boot,1);
MainDecompTypeLabel = cell(Nboots+1,1);
MainDecompTypeLabel{1} = 'MainPCA';
for i = 1:Nboots
    MainDecompTypeLabel{i+1} = ['Bootstrap',num2str(i),'PCA'];
end
ResultsTableHeader = {'BootstrapID','PC1_explained','PC2_explained','PC3_explained','PC4_explained','PC5_explained'};
SourceDataTable = cell2table([MainDecompTypeLabel num2cell([decomp_rand500.explained(1:5)';expl_boot(:,1:5)])],'VariableNames',ResultsTableHeader);
writetable(SourceDataTable,'./SourceDataTables/FigS7a.xlsx')

Nboots = size(expl_boot_tract,1);
TractDecompTypeLabel = cell(Nboots+1,1);
TractDecompTypeLabel{1} = 'TractOnlyPCA';
for i = 1:Nboots
    TractDecompTypeLabel{i+1} = ['BootstrapTractOnly',num2str(i),'PCA'];
end
ResultsTableHeader = {'BootstrapID','PC1_explained','PC2_explained','PC3_explained','PC4_explained','PC5_explained'};
SourceDataTable = cell2table([TractDecompTypeLabel num2cell([tract_decomp.explained(1:5)';expl_boot_tract(:,1:5)])],'VariableNames',ResultsTableHeader);
writetable(SourceDataTable,'./SourceDataTables/FigS7d.xlsx')

SourceData = corr(decomp_rand500.score(:,1),scorePC1_boot)';
SourceDataTable = cell2table([TractDecompTypeLabel(2:end) num2cell(SourceData)],'VariableNames',{'BootstrapID','Correlation_with_PC1score'});
writetable(SourceDataTable,'./SourceDataTables/FigS7b.xlsx')

SourceData = corr(decomp_rand500.coeff(:,1),coeffPC1_boot)';
SourceDataTable = cell2table([TractDecompTypeLabel(2:end) num2cell(SourceData)],'VariableNames',{'BootstrapID','Correlation_with_PC1coeff'});
writetable(SourceDataTable,'./SourceDataTables/FigS7c.xlsx')

SourceData = corr(tract_decomp.score(:,1),scorePC1_boot_tract)';
SourceDataTable = cell2table([TractDecompTypeLabel(2:end) num2cell(SourceData)],'VariableNames',{'BootstrapID','Correlation_with_ConnOnly_PC1score'});
writetable(SourceDataTable,'./SourceDataTables/FigS7e.xlsx')

SourceData = corr(tract_decomp.coeff(:,1),coeffPC1_boot_tract)';
SourceDataTable = cell2table([TractDecompTypeLabel(2:end) num2cell(SourceData)],'VariableNames',{'BootstrapID','Correlation_with_ConnOnly_PC1coeff'});
writetable(SourceDataTable,'./SourceDataTables/FigS7f.xlsx')

%% Plot LOOCV result (Fig S8)
load('loocv_result.mat')

figure
histogram(RMSE)

xlabel('RMSE of PC1 scores')
ylabel({'Number of','leave-one-out iterations'})
set(gca,'FontSize',20)

exportgraphics(gcf,'./figure_outputs/S8.png','Resolution',300)
print('./figure_outputs/S8.svg','-dsvg')

SourceDataTable = array2table(RMSE','VariableNames',{'LOOCV_RMSE'});
writetable(SourceDataTable,'./SourceDataTables/FigS8.xlsx')

%% Plot mouse PCs with coordinates and cortical properties (Figs S9, S10, S19, S20, S21)

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
print('./figure_outputs/S10.svg','-dsvg')
print('./figure_outputs/S10.png','-dpng')

SourceDataTable = cell2table([mThalROIs num2cell([MouseThalNuc_xcoord mouse_pcs_thal(:,1)])],'VariableNames',{'Thalamic_ROI','CCFv3_x_coord','mPC1_score'},'RowNames',mThalROIs);
writetable(SourceDataTable,'./SourceDataTables/FigS10.xlsx')

disp('Mouse PCs with CCFv3 coordinate')

mkdir ./figure_outputs/MouseCoordCorr

PlotMousePCvsCoord(1,{'a','b','c'},'./figure_outputs/MouseCoordCorr/mCoorCorr_');
PlotMousePCvsCoord(2,{'d','e','f'},'./figure_outputs/MouseCoordCorr/mCoorCorr_');
PlotMousePCvsCoord(3,{'g','h','i'},'./figure_outputs/MouseCoordCorr/mCoorCorr_');

SourceDataTable = cell2table([mThalROIs num2cell([MouseThalROICoords mouse_pcs_thal(:,1:3)])],'VariableNames',{'Thalamic_ROI','CCFv3_x_coord','CCFv3_y_coord','CCFv3_z_coord','mPC1_score','mPC2_score','mPC3_score'},'RowNames',mThalROIs);
writetable(SourceDataTable,'./SourceDataTables/FigS9.xlsx')

mouse_CortFeatures = GetAMBAcortexdata();

FigNumber = {'S19','S20','S21'};

ResultsTableHeader = cell(1,size(mouse_CortFeatures.data,2)+1);

% Make Figs S19-S21

for i = 1:3
    PlotMouseCortFeatures(i)
    
    ResultsTableHeader{1} = [FigNumber{i},'_xdata'];
    for j = 1:size(mouse_CortFeatures.data,2)        
        ResultsTableHeader{j+1} = [FigNumber{i},numberToLetter(j),'_ydata'];        
    end    
    
   TableData = [mouse_pcs_cort(:,i) mouse_CortFeatures.data];
   
   SourceDataTable = array2table(TableData,'VariableNames',ResultsTableHeader);
   
   writetable(SourceDataTable,['./SourceDataTables/Fig',FigNumber{i},'.xlsx'])
    
end

%% Compare human and mouse genes for PCs 1-3 (Fig S13)

TractGeneNorm = load('./data/processed/TractGeneNorm_rand500.mat');

GeneNames_human = TractGeneNorm.GeneNames_human;

PhillipsGeneTable = readtable('PhillipsMouseThalGenes.xlsx');
GeneNames_phillips = PhillipsGeneTable.GeneSymbol;
GeneNames_phillips_upper = upper(GeneNames_phillips);

Phillips_pc1_gene = zscore(PhillipsGeneTable.PC1);

[~,idx] = ismember(GeneNames_human,upper(GeneNames_mouse));
AMBA2Human_homolog = idx(idx~=0);
Human2AMBA_homolog = find(idx~=0);
iter = 1;

mkdir ./figure_outputs/MouseVsHumanGene

disp('Human vs Mouse Gene:')

for i = 1:3
    for j = 1:3
    figure
    x = zscore(pcs_gene(Human2AMBA_homolog,i));
    y = zscore(mouse_pcs_gene(AMBA2Human_homolog,j));

    %disp(['Human vs Mouse Gene: Panel ',numberToLetter(iter)])
    DisplayCorr(x,y)

    s = scatterfit(x,y,36,lines(1),[],2,1);
    xlabel(['Human gene PC',num2str(i),' loading'])
    ylabel({'Allen Mouse Brain',['gene PC',num2str(j),' loading']})
    a = annotation('textbox',[0 .89 .05 .13],'String',numberToLetter(iter),'FontSize',32,'EdgeColor','none');
    print(['./figure_outputs/MouseVsHumanGene/MouseVsHumanGene_',num2str(iter),'.png'],'-dpng','-r300')
    print(['./figure_outputs/MouseVsHumanGene/MouseVsHumanGene_',num2str(iter),'.svg'],'-dsvg')
    iter = iter + 1;
    end   
end

SourceDataTable = cell2table([ GeneNames_human(AMBA2Human_homolog) num2cell([zscore(pcs_gene(Human2AMBA_homolog,1:3)) zscore(mouse_pcs_gene(AMBA2Human_homolog,1:3))]) ],'VariableNames',{'Gene','Human_PC1_loading','Human_PC2_loading','Human_PC3_loading','AMBA_PC1_loading','AMBA_PC2_loading','AMBA_PC3_loading'},'RowNames',GeneNames_human(AMBA2Human_homolog));
writetable(SourceDataTable,'./SourceDataTables/FigS13.xlsx')

%% Make neuromap figures (Figs S16, S17, S18)

disp('Getting neuromap figures')

FigNumber = {'S16','S17','S18'};

NeuromapLgd = cell(1,3);

for i = 1:3
    neuromap_corrs = load(['./data/processed/NeuroMapCorrs_PC',num2str(i),'.mat']);
    [SigMapDescrips,SigMapDescripsFull] = PlotSigNeuromaps(neuromap_corrs,['./figure_outputs/SigNeuromaps_PC',num2str(i)],['Cortical region PC',num2str(i),' loading']);
    
    NeuromapLgd{i} = '';
    
    ResultsTableHeader_ = cell(1,length(SigMapDescrips));
    
    for j = 1:length(SigMapDescrips)
        NeuromapLgd{i} = [NeuromapLgd{i},SigMapDescripsFull{j},'. '];        
        ResultsTableHeader_{j} = [FigNumber{i},numberToLetter(j),'_ydata'];        
    end
    
   ResultsTableHeader = [{[FigNumber{i},'_xdata']},ResultsTableHeader_];
    
   TableData = [neuromap_corrs.data neuromap_corrs.neuromap_parc(:,neuromap_corrs.corr_sig)];
   
   SourceDataTable = array2table(TableData,'VariableNames',ResultsTableHeader);
   
   writetable(SourceDataTable,['./SourceDataTables/Fig',FigNumber{i},'.xlsx'])
   
end

%% Makes gene tracjectory plots (Figs, S24, S25, S26)

FigNumber = {'S24','S25','S26'};

for i = 1:3
[NegTrajs,PosTrajs] = plot_gene_trajectories(i);
exportgraphics(gcf,['./figure_outputs/GeneTrajPC',num2str(i),'.png'],'Resolution',300)
print(['./figure_outputs/GeneTrajPC',num2str(i),'.svg'],'-dsvg','-r300')

writetable(NegTrajs,['./SourceDataTables/Fig',FigNumber{i},'a.xlsx'])
writetable(PosTrajs,['./SourceDataTables/Fig',FigNumber{i},'b.xlsx'])

end


filename = './SupplementaryTables/SupplementaryTable10.xlsx';
if exist(filename, 'file')
    delete(filename);
    disp(['File "', filename, '" deleted. Remaking.']);
else
    disp(['Making "', filename,'"']);
end

for i = 1:3
SigDevelopmentalTrajectory = readtable(['./data/processed/gene_trajectories/SigDevelopmentalTrajectoryGenes_PC',num2str(i),'.csv']);
writetable(SigDevelopmentalTrajectory,'./SupplementaryTables/SupplementaryTable10.xlsx','WriteMode','Append','Sheet',['PC',num2str(i),'_','Genes'])
end

close all

%% Make some Source Data Tables

Seeds = (1:length(pcs_thal(:,1)))';

SeedLabel = num2cell(Seeds);

CorticalROIs = (1:length(pcs_cort(:,1)))';

SourceDataTable = array2table([Seeds pcs_thal(:,1)],'VariableNames',{'Thalamic_Seed','PC1_score'});
writetable(SourceDataTable,'./SourceDataTables/Fig2a.xlsx')

SourceDataTable = array2table([Seeds pcs_thal(:,2)],'VariableNames',{'Thalamic_Seed','PC1_score'});
writetable(SourceDataTable,'./SourceDataTables/FigS2a.xlsx')

SourceDataTable = array2table([Seeds pcs_thal(:,3)],'VariableNames',{'Thalamic_Seed','PC1_score'});
writetable(SourceDataTable,'./SourceDataTables/FigS2b.xlsx')


SourceDataTable = array2table([CorticalROIs pcs_cort(:,1)],'VariableNames',{'Cortical_ROI','PC1_loading'});
writetable(SourceDataTable,'./SourceDataTables/Fig3a.xlsx')

SourceDataTable = array2table([CorticalROIs pcs_cort(:,2)],'VariableNames',{'Cortical_ROI','PC1_loading'});
writetable(SourceDataTable,'./SourceDataTables/FigS2c.xlsx')

SourceDataTable = array2table([CorticalROIs pcs_cort(:,3)],'VariableNames',{'Cortical_ROI','PC1_loading'});
writetable(SourceDataTable,'./SourceDataTables/FigS2d.xlsx')

% mThalInds = (1:35)+87;
% mThalROIs = structInfo.acronym(mThalInds);
% mCortROIs = structInfo.acronym(1:38);

SourceDataTable = cell2table([mThalROIs num2cell(mouse_pcs_thal(:,1))],'VariableNames',{'Thalamic_ROI','mPC1_score'},'RowNames',mThalROIs);
writetable(SourceDataTable,'./SourceDataTables/Fig2b.xlsx')

SourceDataTable = cell2table([mThalROIs num2cell(mouse_pcs_thal(:,2))],'VariableNames',{'Thalamic_ROI','mPC1_score'},'RowNames',mThalROIs);
writetable(SourceDataTable,'./SourceDataTables/FigS11a.xlsx')

SourceDataTable = cell2table([mThalROIs num2cell(mouse_pcs_thal(:,3))],'VariableNames',{'Thalamic_ROI','mPC1_score'},'RowNames',mThalROIs);
writetable(SourceDataTable,'./SourceDataTables/FigS11b.xlsx')

SourceDataTable = cell2table([mCortROIs num2cell(mouse_pcs_cort(:,1))],'VariableNames',{'Cortical_ROI','mPC1_loading'},'RowNames',mCortROIs);
writetable(SourceDataTable,'./SourceDataTables/Fig3c.xlsx')

SourceDataTable = cell2table([mCortROIs num2cell(mouse_pcs_cort(:,2))],'VariableNames',{'Cortical_ROI','mPC1_loading'},'RowNames',mCortROIs);
writetable(SourceDataTable,'./SourceDataTables/FigS11c.xlsx')

SourceDataTable = cell2table([mCortROIs num2cell(mouse_pcs_cort(:,3))],'VariableNames',{'Cortical_ROI','mPC1_loading'},'RowNames',mCortROIs);
writetable(SourceDataTable,'./SourceDataTables/FigS11d.xlsx')


%% Make Supplementary Table 1

Seeds = (1:length(pcs_thal(:,1)))';

CorticalROIs = (1:length(pcs_cort(:,1)))';

CorticalROILabel = cell(1,length(CorticalROIs));

for i = 1:length(CorticalROIs)
    CorticalROILabel{i} = ['CorticalROI_',num2str(i)];
end

HumanGeneLabel = cell(1,length(GeneNames_human));

for i = 1:length(GeneNames_human)
    HumanGeneLabel{i} = ['Gene_',GeneNames_human{i}];
end

SuppTable1 = cell2table([num2cell(Seeds) num2cell(pca_input)],'VariableNames',[{'ThalamicSeed'} CorticalROILabel HumanGeneLabel]);

ComponentScoreLabel = cell(1,size(score,2));
ComponentLoadingLabel = cell(1,size(score,2));
ComponentExplained = cell(1,size(score,2));

for i = 1:length(ComponentScoreLabel)
    ComponentScoreLabel{i} = ['PC',num2str(i),'_raw_score'];
    ComponentLoadingLabel{i} = ['PC',num2str(i),'_raw_loading'];
    ComponentExplained{i} = ['PC',num2str(i)];
end

% Returning all the PCs results in a file size which is a bit too big.
% Restrict to the top 500.
TopPCs = 1:500;

SuppTable2_1 = cell2table([num2cell(Seeds) num2cell(score(:,TopPCs))],'VariableNames',[{'ThalamicSeed'} ComponentScoreLabel(TopPCs)]);

SuppTable2_2 = cell2table([ [CorticalROILabel';HumanGeneLabel'] num2cell(coeff(:,TopPCs))],'VariableNames',[{'CorticalROI/Gene'} ComponentScoreLabel(TopPCs)]);

SuppTable2_3 = cell2table([ComponentScoreLabel(TopPCs)' num2cell(explained(TopPCs))],'VariableNames',{'Component','VarianceExplained'});

ComponentScoreLabel = cell(1,size(mouse_score,2));
ComponentLoadingLabel = cell(1,size(mouse_score,2));
ComponentExplained = cell(1,size(mouse_score,2));

for i = 1:length(ComponentScoreLabel)
    ComponentScoreLabel{i} = ['mPC',num2str(i),'_raw_score'];
    ComponentLoadingLabel{i} = ['mPC',num2str(i),'_raw_loading'];
    ComponentExplained{i} = ['mPC',num2str(i)];
end

MouseGeneLabel = cell(1,length(GeneNames_mouse));

for i = 1:length(MouseGeneLabel)
    MouseGeneLabel{i} = ['Gene_',GeneNames_mouse{i}];
end

SuppTable3_1 = cell2table([mThalROIs num2cell(mouse_pca_input)],'VariableNames',[{'ThalamicROI'} mCortROIs' MouseGeneLabel]);

SuppTable4_1 = cell2table([mThalROIs num2cell(mouse_score)],'VariableNames',[{'ThalamicROI'} ComponentScoreLabel]);

SuppTable4_2 = cell2table([ [mCortROIs;MouseGeneLabel'] num2cell(mouse_coeff)],'VariableNames',[{'CorticalROI/Gene'} ComponentLoadingLabel]);

SuppTable4_3 = cell2table([ComponentExplained' num2cell(mouse_explained)],'VariableNames',{'Component','VarianceExplained'});


for i = 1:4
    filename = ['./SupplementaryTables/SupplementaryTable',num2str(i),'.xlsx'];
    if exist(filename, 'file')
        delete(filename);
        disp(['File "', filename, '" deleted. Remaking.']);
    else
        disp(['Making "', filename,'"']);
    end
end

writetable(SuppTable1,'./SupplementaryTables/SupplementaryTable1.xlsx','WriteMode','overwrite','Sheet','HumanConnectivityGeneMatrix')
writetable(SuppTable2_1,'./SupplementaryTables/SupplementaryTable2.xlsx','WriteMode','overwrite','Sheet','HumanPCA_raw_scores')
writetable(SuppTable2_2,'./SupplementaryTables/SupplementaryTable2.xlsx','WriteMode','Append','Sheet','HumanPCA_raw_loadings')
writetable(SuppTable2_3,'./SupplementaryTables/SupplementaryTable2.xlsx','WriteMode','Append','Sheet','HumanPCA_variance_explained')
writetable(SuppTable3_1,'./SupplementaryTables/SupplementaryTable3.xlsx','WriteMode','overwrite','Sheet','MouseConnectivityGeneMatrix')
writetable(SuppTable4_1,'./SupplementaryTables/SupplementaryTable4.xlsx','WriteMode','overwrite','Sheet','MousePCA_raw_scores')
writetable(SuppTable4_2,'./SupplementaryTables/SupplementaryTable4.xlsx','WriteMode','Append','Sheet','MousePCA_raw_loadings')
writetable(SuppTable4_3,'./SupplementaryTables/SupplementaryTable4.xlsx','WriteMode','Append','Sheet','MousePCA_variance_explained')

%% Make supplementary table 3

disp('Make supplementary table 6')


filename = './SupplementaryTables/SupplementaryTable6.xlsx';
if exist(filename, 'file')
    delete(filename);
    disp(['File "', filename, '" deleted. Remaking.']);
else
    disp(['Making "', filename,'"']);
end

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

    writetable(NegTable,'./SupplementaryTables/SupplementaryTable6.xlsx','WriteMode','Append','Sheet',['PC',num2str(PC),'_',NegativeName,'Genes'])
    writetable(PosTable,'./SupplementaryTables/SupplementaryTable6.xlsx','WriteMode','Append','Sheet',['PC',num2str(PC),'_',PositiveName,'Genes'])

end

%%

WebGestalt_FOLDERS = {'PC1Neg','PC1Pos','PC2Neg','PC2Pos','PC3Neg','PC3Pos'};
SheetName = {'PC1_Medial','PC1_Lateral','PC2_Dorsal','PC2_Ventral','PC3_Anterior','PC3_Posterior'};

for i = 11:13
    filename = ['./SupplementaryTables/SupplementaryTable',num2str(i),'.xlsx'];
    if exist(filename, 'file')
        delete(filename);
        disp(['File "', filename, '" deleted. Remaking.']);
    else
        disp(['Making "', filename,'"']);
    end
end

for i = 1:2
[DiseaseEnrich,DiseaseGeneTable,DiseaseEnrichResults] = ReadWebGestalt(['./data/Web/',WebGestalt_FOLDERS{i}]);
writetable(DiseaseEnrich,'./SupplementaryTables/SupplementaryTable11.xlsx','WriteMode','overwrite','Sheet',[SheetName{i},'DiseaseEnrich'])
writetable(DiseaseGeneTable,'./SupplementaryTables/SupplementaryTable11.xlsx','WriteMode','Append','Sheet',[SheetName{i},'DiseaseGenes'])
end

for i = 3:4
[DiseaseEnrich,DiseaseGeneTable,DiseaseEnrichResults] = ReadWebGestalt(['./data/Web/',WebGestalt_FOLDERS{i}]);
writetable(DiseaseEnrich,'./SupplementaryTables/SupplementaryTable12.xlsx','WriteMode','overwrite','Sheet',[SheetName{i},'DiseaseEnrich'])
writetable(DiseaseGeneTable,'./SupplementaryTables/SupplementaryTable12.xlsx','WriteMode','Append','Sheet',[SheetName{i},'DiseaseGenes'])
end

for i = 5:6
[DiseaseEnrich,DiseaseGeneTable,DiseaseEnrichResults] = ReadWebGestalt(['./data/Web/',WebGestalt_FOLDERS{i}]);
writetable(DiseaseEnrich,'./SupplementaryTables/SupplementaryTable13.xlsx','WriteMode','overwrite','Sheet',[SheetName{i},'DiseaseEnrich'])
writetable(DiseaseGeneTable,'./SupplementaryTables/SupplementaryTable13.xlsx','WriteMode','Append','Sheet',[SheetName{i},'DiseaseGenes'])
end

close all