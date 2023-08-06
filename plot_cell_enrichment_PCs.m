function plot_cell_enrichment_PCs(PC,RplOutput)

if nargin < 1
    PC = 1;
end

if nargin < 2
    RplOutput = 0;
end

if PC == 1
    AxisTitle = 'Medial-lateral enrichment ratio';
    NegativeName = 'Medial';
    PositiveName = 'Lateral';
elseif PC == 2
    AxisTitle = 'Dorsal-ventral enrichment ratio';
    NegativeName = 'Dorsal';
    PositiveName = 'Ventral';
elseif PC == 3
    AxisTitle = 'Anterior-posterior enrichment ratio';
    NegativeName = 'Anterior';
    PositiveName = 'Posterior';
end

figure('Position',[0 0 1920 1080])
bar_cmap = brewermap(256,'RdBu');

cell_enrich_table = readtable(['cell_enrichment_PC',num2str(PC),'.csv']);

[isSig,sigLevel,pvals_corr]=BF_FDR(cell_enrich_table.p,.05);

isSigInd = find(isSig);

for i = 1:length(isSigInd)
    idx = isSigInd(i);
disp([cell_enrich_table.loading{idx},' ',cell_enrich_table.class{idx},'  enrichment ratio = ',...
    num2str(cell_enrich_table.enrichment(idx)),' pFDR = ',num2str(pvals_corr(idx)),]);
end

Incell = find(contains(cell_enrich_table.loading,'positive'));

CellName = cell_enrich_table.class(Incell);

PosEnrich = cell_enrich_table.enrichment(Incell);

Incell = find(contains(cell_enrich_table.loading,'negative'));

Negenrich = cell_enrich_table.enrichment(Incell)*-1;



Incell = flipud(find(contains(cell_enrich_table.loading,'positive')));

Cell_Class = cell_enrich_table.class(Incell);

Medial_enrichment = cell_enrich_table.enrichment(Incell+1);

Medial_pFDR = pvals_corr(Incell+1);

Lateral_enrichment = cell_enrich_table.enrichment(Incell);

Lateral_pFDR = pvals_corr(Incell);

EnrichTable = table(Cell_Class,Medial_enrichment,Medial_pFDR,Lateral_enrichment,Lateral_pFDR);

EnrichTable.Properties.VariableNames = {'Cell_class',[NegativeName,'_enrichment'],[NegativeName,'_pFDR'],[PositiveName,'_enrichment'],[PositiveName,'_pFDR']};

axes('Position',[0.1885    0.4787    0.2105    0.4795])

b = barh([PosEnrich Negenrich],'stacked','FaceColor','flat');
 
b(1).FaceColor = bar_cmap(32,:);
b(2).FaceColor = bar_cmap(224,:);

yticks(1:length(CellName))
yticklabels(CellName)

xlabel(AxisTitle)
ylabel('Cell class')
set(gca,'FontSize',20)


% Pos = get(gca,'Position');
% annotation('textbox',[0 Pos(2)+Pos(4)-.1 .1 .1],'String','a','FontSize',32,'EdgeColor','none')
    axis_pos = get(gca,'OuterPosition');
    annotYpos = find_point_on_line(axis_pos(2),axis_pos(2)+axis_pos(4),.9);
    annotation('textbox',[axis_pos(1)-.01 annotYpos axis_pos(3).*.1 axis_pos(2)+axis_pos(4)-annotYpos],'String','a','FontSize',32,'EdgeColor','none')

axes('Position',[0.1885    0.1100    0.2105    0.2389])

neuron_enrich_table = readtable(['neuron_enrichment_PC',num2str(PC),'.csv']);

[isSig,sigLevel,pvals_corr]=BF_FDR(neuron_enrich_table.p,.05);

isSigInd = find(isSig);

Incell = flipud(find(contains(neuron_enrich_table.loading,'positive')));

Neuron_Class = neuron_enrich_table.class(Incell);

Lateral_enrichment = neuron_enrich_table.enrichment(Incell);

Lateral_pFDR = pvals_corr(Incell);

Medial_enrichment = neuron_enrich_table.enrichment(Incell+1);

Medial_pFDR = pvals_corr(Incell+1);

NeuronSubtypeEnrichTable = table(Neuron_Class,Medial_enrichment,Medial_pFDR,Lateral_enrichment,Lateral_pFDR);

NeuronSubtypeEnrichTable.Properties.VariableNames = {'Neuron_class',[NegativeName,'_enrichment'],[NegativeName,'_pFDR'],[PositiveName,'_enrichment'],[PositiveName,'_pFDR']};

Incell = find(contains(neuron_enrich_table.loading,'positive'));

CellName = neuron_enrich_table.class(Incell);

PosEnrich = neuron_enrich_table.enrichment(Incell);

Incell = find(contains(neuron_enrich_table.loading,'negative'));

Negenrich = neuron_enrich_table.enrichment(Incell)*-1;

b = barh([PosEnrich Negenrich],'stacked','FaceColor','flat');
 
b(1).FaceColor = bar_cmap(32,:);
b(2).FaceColor = bar_cmap(224,:);

yticks(1:length(CellName))
yticklabels(CellName)

xlabel(AxisTitle)
ylabel('Neuron class')
set(gca,'FontSize',20)

% Pos = get(gca,'Position');
% annotation('textbox',[0.075 Pos(2)+Pos(4)-.1 .1 .1],'String','b','FontSize',32,'EdgeColor','none')
    axis_pos = get(gca,'OuterPosition');
    annotYpos = find_point_on_line(axis_pos(2),axis_pos(2)+axis_pos(4),.9);
    annotation('textbox',[axis_pos(1)-.01 annotYpos axis_pos(3).*.1 axis_pos(2)+axis_pos(4)-annotYpos],'String','b','FontSize',32,'EdgeColor','none')

subClust_enrich_table = readtable(['subcluster_enrichment_PC',num2str(PC),'.csv']);

ClusterName = {'Habenula','Rora','Gad2/Ahi1'};
NsubPerClust = [4 11 11];

LoadingValence = {'positive','negative'};

subClust_enrich = zeros(sum(NsubPerClust),6);

for i = 1:3

Inclust = find(contains(subClust_enrich_table.subcluster,ClusterName{i}).*contains(subClust_enrich_table.loading,LoadingValence{1}));
Cluster = split(subClust_enrich_table.subcluster(Inclust),'-');

ClusterSub_num = cell2mat(cellfun(@str2num,(Cluster(:,2)),'UniformOutput',false));

[ClusterSub_num_sort,ClusterSub_num_sort_ind] = sort(ClusterSub_num);

InclustInd = Inclust(ClusterSub_num_sort_ind);

PosEnrich = [ClusterSub_num_sort subClust_enrich_table.enrichment(InclustInd) subClust_enrich_table.p(InclustInd)];

Inclust = find(contains(subClust_enrich_table.subcluster,ClusterName{i}).*contains(subClust_enrich_table.loading,LoadingValence{2}));
Cluster = split(subClust_enrich_table.subcluster(Inclust),'-');

ClusterSub_num = cell2mat(cellfun(@str2num,(Cluster(:,2)),'UniformOutput',false));

[ClusterSub_num_sort,ClusterSub_num_sort_ind] = sort(ClusterSub_num);

InclustInd = Inclust(ClusterSub_num_sort_ind);

subClust_enrich_data{i} = [PosEnrich subClust_enrich_table.enrichment(InclustInd)*-1 subClust_enrich_table.p(InclustInd)];

end

subClust_enrich = [ones(NsubPerClust(1),1)*1 subClust_enrich_data{1};ones(NsubPerClust(2),1)*2 subClust_enrich_data{2};ones(NsubPerClust(3),1)*3 subClust_enrich_data{3}];

PosRank = (subClust_enrich(:,3));

NegRank = (subClust_enrich(:,5));

MeanRank = PosRank + NegRank;

[~,MeanRank_order] = sort(MeanRank);



%subplot(2,2,[1 3]) 
%axes('Position',[0.1300    0.1100    0.3347    0.8150])
subcluster_plot = axes('Position',[0.5307    0.1308    0.1693    0.8150]);
 b = barh(subClust_enrich(MeanRank_order,[3 5]),'stacked','FaceColor','flat');
 
 bar_cmap = brewermap(256,'RdBu');
 
 size_cmap = size(bar_cmap,1);
%  neg_bar_cmap = bar_cmap(1:round(size_cmap*.33),:);
%  pos_bar_cmap = flipud(bar_cmap(size_cmap:-1:round(size_cmap*.66),:)); 
 pos_bar_cmap = flipud(bar_cmap(1:128,:));
 neg_bar_cmap = (bar_cmap(size_cmap:-1:129,:)); 
 
 neg_enrich_color = MapData2Colors(subClust_enrich(MeanRank_order,5),neg_bar_cmap);
 pos_enrich_color = MapData2Colors(subClust_enrich(MeanRank_order,3),pos_bar_cmap);
 
 for i = 1:26
     b(1).CData(i,:) = pos_enrich_color(i,:);
     b(2).CData(i,:) = neg_enrich_color(i,:);
 end

ClusterName = {'Habenula','Rora','Gad2/Ahi1'};
NsubPerClust = [4 11 11];
Subcluster_names = cell(26,1);
iter = 1;
for i = 1:3
    for j = 1:NsubPerClust(i)
        Subcluster_names{iter} = [ClusterName{i},'-',num2str(j)];
        iter = iter + 1;
    end
end

yticks(1:26)
yticklabels(Subcluster_names(MeanRank_order))

xlabel(AxisTitle)
ylabel('Neuronal subtypes')

set(gca,'FontSize',20)

Lateral_enrichment = subClust_enrich(flipud(MeanRank_order),3);

Medial_enrichment = subClust_enrich(flipud(MeanRank_order),5);

subclusterPvals = [subClust_enrich(flipud(MeanRank_order),4);subClust_enrich(flipud(MeanRank_order),6)];

[isSig,sigLevel,pvals_corr]=BF_FDR(subclusterPvals,.05);

isSigInd = find(isSig);

Lateral_pFDR = pvals_corr(1:length(Lateral_enrichment));
Medial_pFDR = pvals_corr(length(Lateral_enrichment)+1:end);

Neuronal_subcluster = Subcluster_names(flipud(MeanRank_order));

SubclusterEnrichTable = table(Neuronal_subcluster,Medial_enrichment,Medial_pFDR,Lateral_enrichment,Lateral_pFDR);

SubclusterEnrichTable.Properties.VariableNames = {'Neuronal_subcluster',[NegativeName,'_enrichment'],[NegativeName,'_pFDR'],[PositiveName,'_enrichment'],[PositiveName,'_pFDR']};

% Pos = get(gca,'Position');
% annotation('textbox',[0.45 Pos(2)+Pos(4)-.1 .1 .1],'String','c','FontSize',32,'EdgeColor','none')
    axis_pos = get(gca,'OuterPosition');
    annotYpos = find_point_on_line(axis_pos(2),axis_pos(2)+axis_pos(4),.8);
    annotation('textbox',[axis_pos(1)-.01 annotYpos axis_pos(3).*.1 .17],'String','c','FontSize',32,'EdgeColor','none')

% 
% 
tsnedata = readtable('tsne_subdata.xlsx');

SubClustID_char = tsnedata.cx;
ClusterID_SubID = cell2mat(cellfun(@str2num,(split(SubClustID_char,'-')),'UniformOutput',false));

iter = 1;
subclustid = zeros(length(SubClustID_char),1);
NClusters = unique(ClusterID_SubID(:,1));
for i = 1:length(NClusters)
inclust = ClusterID_SubID(:,1)==NClusters(i);
u = unique(ClusterID_SubID(inclust,2))';
for j = u
   insub_clust = ClusterID_SubID(:,2)==j & inclust;
   subclustid(insub_clust) = iter;
   iter = iter + 1;
end
end

% axes('Position',[TSNE_plots_xpos    .55  (1-TSNE_plots_xpos)*.95 .35])

TSNE_plots_xpos_neg = subcluster_plot.Position(1)+subcluster_plot.Position(3)+.01;

TSNE_plots_width = (.99-TSNE_plots_xpos_neg)/2;

TSNE_plots_total_height = subcluster_plot.Position(4);

TSNE_plots_height = (TSNE_plots_total_height)/3;

TSNE_plots_xpos_pos = TSNE_plots_xpos_neg+TSNE_plots_width;

TSNE_plots_height_start = subcluster_plot.Position(2);

neg_enrich_limts = [min(subClust_enrich(:,5)) max(subClust_enrich(:,5))];

pos_enrich_limts = [min(subClust_enrich(:,3)) max(subClust_enrich(:,3))];

run = 1;
if run

for i = 1:3

axes('Position',[TSNE_plots_xpos_pos   TSNE_plots_height_start+(TSNE_plots_height*(i-1))    TSNE_plots_width    TSNE_plots_height])
    
insubclust = ClusterID_SubID(:,1)==i;

insubclust_id = subclustid(insubclust);

Enrich_val = changem(insubclust_id,subClust_enrich(:,3),1:26);

Enrich_val_color = MapData2Colors(Enrich_val,pos_bar_cmap,pos_enrich_limts);

%Enrich_val_color(Enrich_val==0,:) = 1;

scatter(tsnedata.V1(insubclust),tsnedata.V2(insubclust),30,Enrich_val_color,'filled','MarkerEdgeColor','k')

axis tight
axis off

ylimits = ylim;
new_ymin = find_point_on_line(ylimits(2),ylimits(1),1.15);
new_ymax = find_point_on_line(ylimits(2),ylimits(1),-0.05);
text_y = find_point_on_line(ylimits(2),ylimits(1),1.075);
ylim([new_ymin new_ymax])
xlimits = xlim;

new_xmin = find_point_on_line(xlimits(2),xlimits(1),1.05);
new_xmax = find_point_on_line(xlimits(2),xlimits(1),-0.05);
xlim([new_xmin new_xmax])
t = text(mean([new_xmin new_xmax]),text_y,ClusterName{i},'FontSize',20,'HorizontalAlignment','center');
end

for i = 1:3
axes('Position',[TSNE_plots_xpos_neg   TSNE_plots_height_start+(TSNE_plots_height*(i-1))    TSNE_plots_width    TSNE_plots_height])
    
insubclust = ClusterID_SubID(:,1)==i;

insubclust_id = subclustid(insubclust);

Enrich_val = changem(insubclust_id,subClust_enrich(:,5),1:26);

Enrich_val_color = MapData2Colors(Enrich_val,neg_bar_cmap,neg_enrich_limts);

%Enrich_val_color(Enrich_val==0,:) = 1;

scatter(tsnedata.V1(insubclust),tsnedata.V2(insubclust),30,Enrich_val_color,'filled','MarkerEdgeColor','k')

axis tight
axis off

ylimits = ylim;
new_ymin = find_point_on_line(ylimits(2),ylimits(1),1.15);
new_ymax = find_point_on_line(ylimits(2),ylimits(1),-0.05);
text_y = find_point_on_line(ylimits(2),ylimits(1),1.075);
ylim([new_ymin new_ymax])
xlimits = xlim;

new_xmin = find_point_on_line(xlimits(2),xlimits(1),1.05);
new_xmax = find_point_on_line(xlimits(2),xlimits(1),-0.05);
xlim([new_xmin new_xmax])
t = text(mean([new_xmin new_xmax]),text_y,ClusterName{i},'FontSize',20,'HorizontalAlignment','center');
end

    axis_pos = get(gca,'OuterPosition');
    annotYpos = find_point_on_line(axis_pos(2),axis_pos(2)+axis_pos(4),.9);
    annotation('textbox',[axis_pos(1)+.015 annotYpos axis_pos(3).*.1 axis_pos(2)+axis_pos(4)-annotYpos],'String','d','FontSize',32,'EdgeColor','none')

end
exportgraphics(gcf,['./figure_outputs/cell_enrich_PC',num2str(PC),'.png'],'Resolution',300)

if RplOutput
 
delete(['CellEnrichmentCombined_PC',num2str(PC),'.xlsx'])

writetable(EnrichTable,['CellEnrichmentCombined_PC',num2str(PC),'.xlsx'],'WriteMode','overwrite','Sheet','CellClassEnrichment')

writetable(NeuronSubtypeEnrichTable,['CellEnrichmentCombined_PC',num2str(PC),'.xlsx'],'WriteMode','Append','Sheet','NeuronSubtypeEnrichment')

writetable(SubclusterEnrichTable,['CellEnrichmentCombined_PC',num2str(PC),'.xlsx'],'WriteMode','Append','Sheet','NeuronSubclusterEnrichment')

NeuronSubtypes={'Habenula_Tac2','Rora','Gad2-Ahi1'};

for i = 1:3
    intable=readtable(['NegativeNeuron_',NeuronSubtypes{i},'Genes_PC',num2str(PC),'.csv']);
    inGenes{i}=intable.Gene;
end
for i = 1:3
    intable=readtable(['PositiveNeuron_',NeuronSubtypes{i},'Genes_PC',num2str(PC),'.csv']);
    inGenes{i+3}=intable.Gene;
end

for i = 1:6
    Ngenes(i) = length(inGenes{i});
end

MaxGenes = max(Ngenes);

for i = 1:6
    if Ngenes(i) < MaxGenes
    inGenes{i}(Ngenes(i)+1:MaxGenes) = {''};
    end
end

for i = 1:6
    SubtybeGenesPadded(:,i) = inGenes{i};
end

for i = 1:3
NeuronSubtypeGeneVarName{i} = [NegativeName,'_',NeuronSubtypes{i}];
NeuronSubtypeGeneVarName{i+3} = [PositiveName,'_',NeuronSubtypes{i}];
end

NeuronClassGenes = cell2table(SubtybeGenesPadded,'VariableNames',NeuronSubtypeGeneVarName);

writetable(NeuronClassGenes,['CellEnrichmentCombined_PC',num2str(PC),'.xlsx'],'WriteMode','Append','Sheet','NeuronClassGenes')

end