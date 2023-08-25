function plot_cell_enrichment_PCs(PC,TableOutput)

% Check if PC argument is provided, otherwise default to PC 1
if nargin < 1
    PC = 1;
end

% Initialize the MakeTable variable to control table creation
MakeTable = 1;

% Check if TableOutput argument is provided, otherwise skip table creation
if nargin < 2
    MakeTable = 0;
end

% Determine the axis title and names based on the selected principal component (PC)
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
elseif PC >= 4
    AxisTitle = 'Negative-positive enrichment ratio';
    NegativeName = 'Negative';
    PositiveName = 'Positive';
end

% Create a new figure with specific dimensions
figure('Position', [0 0 1920 1080])

% Create a color map for the bar plot
bar_cmap = brewermap(256, 'RdBu');

% Read the cell enrichment table for the selected principal component (PC)
cell_enrich_table = readtable(['cell_enrichment_PC', num2str(PC), '.csv']);

% Perform statistical analysis for significance
[isSig, sigLevel, pvals_corr] = BF_FDR(cell_enrich_table.p, .05);

% Identify indices of significant enrichments
isSigInd = find(isSig);

% Display information about significant enrichments
for i = 1:length(isSigInd)
    idx = isSigInd(i);
    disp([cell_enrich_table.loading{idx}, ' ', cell_enrich_table.class{idx},...
        '  enrichment ratio = ', num2str(cell_enrich_table.enrichment(idx)),...
        ' pFDR = ', num2str(pvals_corr(idx)),]);
end

% Identify indices of positive enrichments
Incell = find(contains(cell_enrich_table.loading, 'positive'));
CellName = cell_enrich_table.class(Incell);
PosEnrich = cell_enrich_table.enrichment(Incell);

% Identify indices of negative enrichments
Incell = find(contains(cell_enrich_table.loading, 'negative'));
Negenrich = cell_enrich_table.enrichment(Incell) * -1;

%%

% Identify indices of positive enrichments and flip their order
Incell = flipud(find(contains(cell_enrich_table.loading, 'positive')));

% Extract relevant information for positive enrichments
Cell_Class = cell_enrich_table.class(Incell);
Medial_enrichment = cell_enrich_table.enrichment(Incell + 1);
Medial_pFDR = pvals_corr(Incell + 1);
Lateral_enrichment = cell_enrich_table.enrichment(Incell);
Lateral_pFDR = pvals_corr(Incell);

% Create a table to store enriched information
EnrichTable = table(Cell_Class, Medial_enrichment, Medial_pFDR, Lateral_enrichment, Lateral_pFDR);

% Set custom variable names for the enriched table
EnrichTable.Properties.VariableNames = {'Cell_class', ...
    [NegativeName, '_enrichment'], [NegativeName, '_pFDR'], ...
    [PositiveName, '_enrichment'], [PositiveName, '_pFDR']};

% Define the position of the axes for the enrichment bar plots

CellClassAxesPos = [0.1641    0.4352    0.2146    0.537];

NeuronClassAxesPos = [0.1641    0.0759    0.2146    0.2547];

NeuronSubtypeAxesPos = [0.5141    0.0759    0.1948    CellClassAxesPos(2)+CellClassAxesPos(4)-0.0759];

% Create a horizontal bar plot for cell class enrichment
axes('Position', CellClassAxesPos)

% Create stacked bar plot with positive and negative enrichments
b = barh([PosEnrich, Negenrich], 'stacked', 'FaceColor', 'flat');

% Set the face colors for the stacked bars
b(1).FaceColor = bar_cmap(32, :);
b(2).FaceColor = bar_cmap(224, :);

% Set y-axis ticks and labels
yticks(1:length(CellName))
yticklabels(CellName)

% Set x-axis and y-axis labels
xlabel(AxisTitle)
ylabel('Cell class')

% Set font size for the current axis
set(gca, 'FontSize', 20)

%%
% Retrieve the current axis' outer position
axis_pos = get(gca, 'OuterPosition');

% Calculate the vertical position for the annotation 'a'
% along the specified line fraction
annotYpos = find_point_on_line(axis_pos(2), axis_pos(2) + axis_pos(4), 0.835);

% Define the position for annotation 'a'
aPanelLoc = [axis_pos(1), annotYpos, axis_pos(3) * 0.1, 0.1];

% Create a text annotation 'a' with specified position and appearance
annotation('textbox', aPanelLoc, 'String', 'a', 'FontSize', 32, 'EdgeColor', 'none');

% Create new axes for neuron class enrichment bar plot using defined position
axes('Position', NeuronClassAxesPos)

% Read neuron enrichment table from CSV file based on the PC value
neuron_enrich_table = readtable(['neuron_enrichment_PC', num2str(PC), '.csv']);

% Perform Benjamini-Hochberg False Discovery Rate correction
[isSig, sigLevel, pvals_corr] = BF_FDR(neuron_enrich_table.p, 0.05);

% Find the indices of significant enrichments
isSigInd = find(isSig);

% Identify positive enrichments and flip their order
Incell = flipud(find(contains(neuron_enrich_table.loading, 'positive')));

% Extract relevant information for positive enrichments
Neuron_Class = neuron_enrich_table.class(Incell);
Lateral_enrichment = neuron_enrich_table.enrichment(Incell);
Lateral_pFDR = pvals_corr(Incell);
Medial_enrichment = neuron_enrich_table.enrichment(Incell + 1);
Medial_pFDR = pvals_corr(Incell + 1);

% Create a table to store neuron subtype enrichment information
NeuronSubtypeEnrichTable = table(Neuron_Class, Medial_enrichment, Medial_pFDR, Lateral_enrichment, Lateral_pFDR);

% Set custom variable names for the neuron subtype enrichment table
NeuronSubtypeEnrichTable.Properties.VariableNames = {'Neuron_class', ...
    [NegativeName, '_enrichment'], [NegativeName, '_pFDR'], ...
    [PositiveName, '_enrichment'], [PositiveName, '_pFDR']};

% Identify positive enrichments for neuron class
Incell = find(contains(neuron_enrich_table.loading, 'positive'));

% Extract relevant information for positive enrichments
CellName = neuron_enrich_table.class(Incell);
PosEnrich = neuron_enrich_table.enrichment(Incell);

% Identify negative enrichments for neuron class
Incell = find(contains(neuron_enrich_table.loading, 'negative'));

% Calculate negative enrichments as negative values
Negenrich = neuron_enrich_table.enrichment(Incell) * -1;

% Create a horizontal bar plot for neuron class enrichment
b = barh([PosEnrich, Negenrich], 'stacked', 'FaceColor', 'flat');

% Set the face colors for the stacked bars
b(1).FaceColor = bar_cmap(32, :);
b(2).FaceColor = bar_cmap(224, :);

% Set y-axis ticks and labels
yticks(1:length(CellName))
yticklabels(CellName)

% Set x-axis and y-axis labels
xlabel(AxisTitle)
ylabel('Neuron class')

% Set font size for the current axis
set(gca, 'FontSize', 20)

%%

% Retrieve the outer position of the current axis
axis_pos = get(gca, 'OuterPosition');

% Calculate the vertical position for the annotation 'b'
annotYpos = find_point_on_line(axis_pos(2), axis_pos(2) + axis_pos(4), 0.9);

% Define the position for annotation 'b'
annotation_x = axis_pos(1);
annotation_width = axis_pos(3) * 0.1;
annotation_y = annotYpos;
annotation_height = axis_pos(2) + axis_pos(4) - annotYpos;
annotation_position = [annotation_x, annotation_y, annotation_width, annotation_height];

% Create a text annotation 'b' with specified position and appearance
annotation('textbox', annotation_position, 'String', 'b', 'FontSize', 32, 'EdgeColor', 'none');

% Read subcluster enrichment data from a CSV file based on the PC value
subClust_enrich_table = readtable(['subcluster_enrichment_PC', num2str(PC), '.csv']);

% Define the names of the subclusters and the number of subclusters per cluster
ClusterName = {'Habenula', 'Rora', 'Gad2/Ahi1'};
NsubPerClust = [4, 11, 11];

% Define loading valence categories
LoadingValence = {'positive', 'negative'};

% Initialize an empty matrix to store subcluster enrichment data
subClust_enrich = zeros(sum(NsubPerClust), 6);

% Loop through each cluster
for i = 1:3
    % Find indices of subclusters with specific loading valence for the current cluster
    Inclust = find(contains(subClust_enrich_table.subcluster, ClusterName{i}) .* contains(subClust_enrich_table.loading, LoadingValence{1}));
    Cluster = split(subClust_enrich_table.subcluster(Inclust), '-');

    % Extract and sort subcluster numbers
    ClusterSub_num = cell2mat(cellfun(@str2num, (Cluster(:, 2)), 'UniformOutput', false));
    [ClusterSub_num_sort, ClusterSub_num_sort_ind] = sort(ClusterSub_num);
    InclustInd = Inclust(ClusterSub_num_sort_ind);

    % Extract enrichment data for positive valence and create a matrix
    PosEnrich = [ClusterSub_num_sort, subClust_enrich_table.enrichment(InclustInd), subClust_enrich_table.p(InclustInd)];

    % Find indices of subclusters with the opposite loading valence for the current cluster
    Inclust = find(contains(subClust_enrich_table.subcluster, ClusterName{i}) .* contains(subClust_enrich_table.loading, LoadingValence{2}));
    Cluster = split(subClust_enrich_table.subcluster(Inclust), '-');

    % Extract and sort subcluster numbers
    ClusterSub_num = cell2mat(cellfun(@str2num, (Cluster(:, 2)), 'UniformOutput', false));
    [ClusterSub_num_sort, ClusterSub_num_sort_ind] = sort(ClusterSub_num);
    InclustInd = Inclust(ClusterSub_num_sort_ind);

    % Extract enrichment data for negative valence and create a matrix
    subClust_enrich_data{i} = [PosEnrich, subClust_enrich_table.enrichment(InclustInd) * -1, subClust_enrich_table.p(InclustInd)];
end

% Combine the subcluster enrichment data into a single matrix
subClust_enrich = [ones(NsubPerClust(1), 1) * 1, subClust_enrich_data{1};
                   ones(NsubPerClust(2), 1) * 2, subClust_enrich_data{2};
                   ones(NsubPerClust(3), 1) * 3, subClust_enrich_data{3}];

% Calculate positive and negative rank sums
PosRank = (subClust_enrich(:, 3));
NegRank = (subClust_enrich(:, 5));

% Calculate the mean rank by summing positive and negative ranks
MeanRank = PosRank + NegRank;

% Sort the subclusters by mean rank
[~, MeanRank_order] = sort(MeanRank);

%%
% Create a new axes for the subcluster plot
subcluster_plot = axes('Position', NeuronSubtypeAxesPos);

% Create a horizontal stacked bar plot for subcluster enrichment
b = barh(subClust_enrich(MeanRank_order, [3, 5]), 'stacked', 'FaceColor', 'flat');

% Generate a colormap for the bar colors
bar_cmap = brewermap(256, 'RdBu');
size_cmap = size(bar_cmap, 1);
pos_bar_cmap = flipud(bar_cmap(1:128, :));
neg_bar_cmap = (bar_cmap(size_cmap:-1:129, :)); 

% Map enrichment values to colors for positive and negative enrichments
neg_enrich_color = MapData2Colors(subClust_enrich(MeanRank_order, 5), neg_bar_cmap);
pos_enrich_color = MapData2Colors(subClust_enrich(MeanRank_order, 3), pos_bar_cmap);

% Apply the calculated colors to the bar plot
for i = 1:26
    b(1).CData(i, :) = pos_enrich_color(i, :);
    b(2).CData(i, :) = neg_enrich_color(i, :);
end

% Define subcluster names based on cluster names and subcluster numbers
ClusterName = {'Habenula', 'Rora', 'Gad2/Ahi1'};
NsubPerClust = [4, 11, 11];
Subcluster_names = cell(26, 1);
iter = 1;
for i = 1:3
    for j = 1:NsubPerClust(i)
        Subcluster_names{iter} = [ClusterName{i}, '-', num2str(j)];
        iter = iter + 1;
    end
end

% Customize y-axis ticks and labels
yticks(1:26)
yticklabels(Subcluster_names(MeanRank_order))

% Set labels and font size for the plot
xlabel(AxisTitle)
ylabel('Neuronal subtypes')
set(gca, 'FontSize', 20)

% Extract enrichment data for lateral and medial subclusters
Lateral_enrichment = subClust_enrich(flipud(MeanRank_order), 3);
Medial_enrichment = subClust_enrich(flipud(MeanRank_order), 5);

% Combine p-values for lateral and medial subclusters
subclusterPvals = [subClust_enrich(flipud(MeanRank_order), 4); subClust_enrich(flipud(MeanRank_order), 6)];

% Perform Benjamini-Hochberg FDR correction
[isSig, sigLevel, pvals_corr] = BF_FDR(subclusterPvals, 0.05);

% Identify significant indices
isSigInd = find(isSig);

% Separate corrected p-values for lateral and medial subclusters
Lateral_pFDR = pvals_corr(1:length(Lateral_enrichment));
Medial_pFDR = pvals_corr(length(Lateral_enrichment)+1:end);

% Extract neuronal subcluster names
Neuronal_subcluster = Subcluster_names(flipud(MeanRank_order));

% Create a table for subcluster enrichment data
SubclusterEnrichTable = table(Neuronal_subcluster, Medial_enrichment, Medial_pFDR, Lateral_enrichment, Lateral_pFDR);

% Customize variable names for the table
SubclusterEnrichTable.Properties.VariableNames = {'Neuronal_subcluster', [NegativeName, '_enrichment'], [NegativeName, '_pFDR'], [PositiveName, '_enrichment'], [PositiveName, '_pFDR']};

% Retrieve the outer position of the current axis
axis_pos = get(gca, 'OuterPosition');

% Calculate the vertical position for the annotation 'c'
annotYpos = find_point_on_line(axis_pos(2), axis_pos(2) + axis_pos(4), 0.9);

% Create annotation 'c' with specified position and appearance
annotation('textbox', [axis_pos(1), aPanelLoc(2), axis_pos(3) * 0.1, aPanelLoc(4)], 'String', 'c', 'FontSize', 32, 'EdgeColor', 'none');

% Read t-SNE data from an Excel file
tsnedata = readtable('tsne_subdata.xlsx');

% Extract cluster and subcluster IDs from the data
SubClustID_char = tsnedata.cx;
ClusterID_SubID = cell2mat(cellfun(@str2num, (split(SubClustID_char, '-')), 'UniformOutput', false));

% Initialize variables for subcluster identification
iter = 1;
subclustid = zeros(length(SubClustID_char), 1);
NClusters = unique(ClusterID_SubID(:, 1));

% Assign unique IDs to subclusters based on cluster and subcluster IDs
for i = 1:length(NClusters)
    inclust = ClusterID_SubID(:, 1) == NClusters(i);
    u = unique(ClusterID_SubID(inclust, 2))';
    for j = u
        insub_clust = ClusterID_SubID(:, 2) == j & inclust;
        subclustid(insub_clust) = iter;
        iter = iter + 1;
    end
end

% Define positions for the t-SNE plot axes
TSNE_plots_xpos_neg = subcluster_plot.Position(1) + subcluster_plot.Position(3) + 0.01;
TSNE_plots_width = (0.9999 - TSNE_plots_xpos_neg) / 2;
TSNE_plots_total_height = subcluster_plot.Position(4);
TSNE_plots_height = TSNE_plots_total_height / 3;
TSNE_plots_xpos_pos = TSNE_plots_xpos_neg + TSNE_plots_width;
TSNE_plots_height_start = subcluster_plot.Position(2);

% Define enrichment limits for negative and positive enrichments
neg_enrich_limts = [min(subClust_enrich(:, 5)), max(subClust_enrich(:, 5))];
pos_enrich_limts = [min(subClust_enrich(:, 3)), max(subClust_enrich(:, 3))];

% Check if the "run" flag is true
run = 1;
if run
    
    % Loop through each cluster
    for i = 1:3
        
        % Create axes for t-SNE plot on positive enrichment side
        axes('Position', [TSNE_plots_xpos_pos, TSNE_plots_height_start + (TSNE_plots_height * (i - 1)), TSNE_plots_width, TSNE_plots_height])
        
        % Identify subclusters within the current cluster
        insubclust = ClusterID_SubID(:, 1) == i;
        insubclust_id = subclustid(insubclust);
        
        % Map enrichment values to colors
        Enrich_val = changem(insubclust_id, subClust_enrich(:, 3), 1:26);
        Enrich_val_color = MapData2Colors(Enrich_val, pos_bar_cmap, pos_enrich_limts);
        
        % Create a scatter plot with the mapped colors
        scatter(tsnedata.V1(insubclust), tsnedata.V2(insubclust), 30, Enrich_val_color, 'filled', 'MarkerEdgeColor', 'k')
        
        % Customize plot appearance
        axis tight
        axis off
        ylimits = ylim;
        new_ymin = find_point_on_line(ylimits(2), ylimits(1), 1.15);
        new_ymax = find_point_on_line(ylimits(2), ylimits(1), -0.05);
        text_y = find_point_on_line(ylimits(2), ylimits(1), 1.075);
        ylim([new_ymin new_ymax])
        xlimits = xlim;
        new_xmin = find_point_on_line(xlimits(2), xlimits(1), 1.05);
        new_xmax = find_point_on_line(xlimits(2), xlimits(1), -0.05);
        xlim([new_xmin new_xmax])
        
        % Add text label indicating the cluster name
        t = text(mean([new_xmin new_xmax]), text_y, ClusterName{i}, 'FontSize', 20, 'HorizontalAlignment', 'center');
    end
    
    % Loop through each cluster for negative enrichment side
    for i = 1:3
        
        % Create axes for t-SNE plot on negative enrichment side
        axes('Position', [TSNE_plots_xpos_neg, TSNE_plots_height_start + (TSNE_plots_height * (i - 1)), TSNE_plots_width, TSNE_plots_height])
        
        % Identify subclusters within the current cluster
        insubclust = ClusterID_SubID(:, 1) == i;
        insubclust_id = subclustid(insubclust);
        
        % Map enrichment values to colors
        Enrich_val = changem(insubclust_id, subClust_enrich(:, 5), 1:26);
        Enrich_val_color = MapData2Colors(Enrich_val, neg_bar_cmap, neg_enrich_limts);
        
        % Create a scatter plot with the mapped colors
        scatter(tsnedata.V1(insubclust), tsnedata.V2(insubclust), 30, Enrich_val_color, 'filled', 'MarkerEdgeColor', 'k')
        
        % Customize plot appearance
        axis tight
        axis off
        ylimits = ylim;
        new_ymin = find_point_on_line(ylimits(2), ylimits(1), 1.15);
        new_ymax = find_point_on_line(ylimits(2), ylimits(1), -0.05);
        text_y = find_point_on_line(ylimits(2), ylimits(1), 1.075);
        ylim([new_ymin new_ymax])
        xlimits = xlim;
        new_xmin = find_point_on_line(xlimits(2), xlimits(1), 1.05);
        new_xmax = find_point_on_line(xlimits(2), xlimits(1), -0.05);
        xlim([new_xmin new_xmax])
        
        % Add text label indicating the cluster name
        t = text(mean([new_xmin new_xmax]), text_y, ClusterName{i}, 'FontSize', 20, 'HorizontalAlignment', 'center');
    end
    
    % Calculate position for annotation 'd'
    axis_pos = get(gca, 'OuterPosition');
    annotation('textbox', [axis_pos(1) + 0.015, aPanelLoc(2), axis_pos(3) * 0.1, aPanelLoc(4)], 'String', 'd', 'FontSize', 32, 'EdgeColor', 'none');
end

%% Make tables

% Check if MakeTable flag is true
if MakeTable
    
    % Remove the output file if it exists because this causes issues when
    % appending sheets
    filename = TableOutput;
    if exist(filename, 'file')
        delete(filename);
        disp(['File "', filename, '" deleted.']);
    else
        disp(['File "', filename, '" does not already exist']);
    end
    
    % Write EnrichTable to the specified Excel file, overwriting if necessary
    writetable(EnrichTable, TableOutput, 'WriteMode', 'overwrite', 'Sheet', 'CellClassEnrichment')
    
    % Append NeuronSubtypeEnrichTable to the specified Excel file
    writetable(NeuronSubtypeEnrichTable, TableOutput, 'WriteMode', 'Append', 'Sheet', 'NeuronSubtypeEnrichment')
    
    % Append SubclusterEnrichTable to the specified Excel file
    writetable(SubclusterEnrichTable, TableOutput, 'WriteMode', 'Append', 'Sheet', 'NeuronSubclusterEnrichment')
    
    % Define NeuronSubtypes
    NeuronSubtypes = {'Habenula_Tac2', 'Rora', 'Gad2-Ahi1'};
    
    % Loop to read genes from NegativeNeuron CSV files for each NeuronSubtype
    for i = 1:3
        intable = readtable(['NegativeNeuron_', NeuronSubtypes{i}, 'Genes_PC', num2str(PC), '.csv']);
        inGenes{i} = intable.Gene;
    end
    
    % Loop to read genes from PositiveNeuron CSV files for each NeuronSubtype
    for i = 1:3
        intable = readtable(['PositiveNeuron_', NeuronSubtypes{i}, 'Genes_PC', num2str(PC), '.csv']);
        inGenes{i + 3} = intable.Gene;
    end
    
    % Loop to calculate the number of genes for each subtype
    for i = 1:6
        Ngenes(i) = length(inGenes{i});
    end
    
    % Find the maximum number of genes among all subtypes
    MaxGenes = max(Ngenes);
    
    % Pad gene lists to have the same length for writing to Excel
    for i = 1:6
        if Ngenes(i) < MaxGenes
            inGenes{i}(Ngenes(i) + 1:MaxGenes) = {''};
        end
    end
    
    % Combine padded gene lists into a table
    for i = 1:6
        SubtybeGenesPadded(:, i) = inGenes{i};
    end
    
    % Define variable names for the NeuronClassGenes table
    for i = 1:3
        NeuronSubtypeGeneVarName{i} = [NegativeName, '_', NeuronSubtypes{i}];
        NeuronSubtypeGeneVarName{i + 3} = [PositiveName, '_', NeuronSubtypes{i}];
    end
    
    % Create a table with padded gene lists and specified variable names
    NeuronClassGenes = cell2table(SubtybeGenesPadded, 'VariableNames', NeuronSubtypeGeneVarName);
    
    % Append NeuronClassGenes table to the specified Excel file
    writetable(NeuronClassGenes, TableOutput, 'WriteMode', 'Append', 'Sheet', 'NeuronClassGenes')
end
