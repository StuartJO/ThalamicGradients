function [NegTrajs, PosTrajs] = plot_gene_trajectories(PC)

% This function plots gene expression trajectories based on principal components.

% Inputs:
%   PC: Principal component index for which gene expression trajectories are plotted.
%       Default value is 1 if not provided.

% Outputs:
%   NegTrajs: A table containing negative gene expression trajectories. Each column
%             corresponds to a gene, and the first column contains age values.
%   PosTrajs: A table containing positive gene expression trajectories. Each column
%             corresponds to a gene, and the first column contains age values.

% Check if PC input is provided, otherwise default to PC1
if nargin < 1
    PC = 1;
end

% Create color maps for negative and positive trajectories
bar_cmap = brewermap(256, 'RdBu');
size_cmap = size(bar_cmap, 1);
pos_bar_cmap = flipud(bar_cmap(1:128, :));
neg_bar_cmap = (bar_cmap(129:size_cmap, :));
RedColor = bar_cmap(32, :);
BlueColor = bar_cmap(224, :);

% Create a figure for plotting
figure('Position', [68 192 1317 641])

% Loop for negative and positive trajectories
for i = 1:2
    subplot(1, 2, i)
    if i == 1
        % Load negative trajectories data
        traj = readtable(['negative_trajectories_PC', num2str(PC), '.csv']);
        plotcmap = neg_bar_cmap;
        Plotlabel = 'a';
    else
        % Load positive trajectories data
        traj = readtable(['positive_trajectories_PC', num2str(PC), '.csv']);
        plotcmap = pos_bar_cmap;
        Plotlabel = 'b';
    end

    % Get unique gene symbols
    GenestrajID = unique(traj.symbol);

    % Define the range for enrichment values
    enrichlimits = [min(traj.fit), max(traj.fit)];

    % Calculate the number of unique ages
    Nages = length(unique(traj.age));

    % Initialize a matrix to store fitted values
    FITS = zeros(Nages, length(GenestrajID));

    hold on
    for j = 1:length(GenestrajID)
        IND = strcmp(traj.symbol, GenestrajID{j});
        x = traj.age(IND);
        y = traj.fit(IND);

        % Map fitted values to colormap
        plotcolor = MapData2Colors(max(y), plotcmap, enrichlimits);

        % Plot gene expression trajectory
        pos_traj_plot(j) = plot(x, y, 'Color', [plotcolor, 1]);

        FITS(:, j) = y;

    end

    AGE = x;

    if i == 1
        % Create the NegTrajs table
        NegTrajs = array2table([AGE, FITS], 'VariableNames', [{'Age_in_days'}, GenestrajID']);
    else
        % Create the PosTrajs table
        PosTrajs = array2table([AGE, FITS], 'VariableNames', [{'Age_in_days'}, GenestrajID']);
    end

    % Set the x-axis to log scale
    set(gca, 'XScale', 'log')

    % Adjust y-axis limits
    ylimits = ylim;
    ylim([ylimits(1), find_point_on_line(ylimits(1), ylimits(2), 1.1)])
    ylimits = ylim;

    % Add a vertical line at postnatal age
    plot([280, 280], ylimits, 'Color', 'k')

    % Add labels for prenatal and postnatal
    t1 = text(280, find_point_on_line(ylimits(1), ylimits(2), .95), '  postnatal', 'HorizontalAlignment', 'left', 'FontSize', 20);
    t2 = text(280, find_point_on_line(ylimits(1), ylimits(2), .95), 'prenatal  ', 'HorizontalAlignment', 'right', 'FontSize', 20);

    % Add axis labels and title
    xlabel('Age in days (log)')
    ylabel('Relative gene expression')
    title(['PC', num2str(PC)])
    set(gca, 'FontSize', 20)

    % Add subplot annotation
    axis_pos = get(gca, 'OuterPosition');
    annotYpos = find_point_on_line(axis_pos(2), axis_pos(2) + axis_pos(4), .9);
    annotation('textbox', [axis_pos(1) + .015, annotYpos, axis_pos(3) * .1, axis_pos(2) + axis_pos(4) - annotYpos], 'String', Plotlabel, 'FontSize', 24, 'EdgeColor', 'none')
end
