function FigCaption = PlotMousePCvsCoord(grad, annotlabels, SavePreFix)
% Generates scatter plots of mouse PC scores against CCFv3 coordinates and displays correlation statistics.
%
% Inputs:
%   grad - Principal component (PC) gradient index
%   annotlabels - Labels for the spatial axes (e.g., {'X', 'Y', 'Z'})
%   SavePreFix - Prefix for saving the generated figures
%
% Outputs:
%   FigCaption - Cell array containing figure captions
%
% Usage example:
%   FigCaption = PlotMousePCvsCoord(grad, {'X', 'Y', 'Z'}, 'Mouse_PC_vs_Coord');
%
%   This function generates scatter plots of mouse PC scores against CCFv3 coordinates
%   along different spatial axes and displays correlation statistics.
% Load necessary data
SPATIAL_AXIS = {'x', 'y', 'z'};
SPATIAL_DIR = {'Medial-lateral', 'Anterior-posterior', 'Ventral-dorsal'};
cmap = turbo(256);

load('./data/processed/mouse_decomp.mat', 'mouse_pcs_thal')
load('MouseThalROICoords.mat', 'MouseThalROICoords')

FigCaption = cell(3, 1);

mpc_thal = zscore(mouse_pcs_thal(:, grad));
for XYZ = 1:3
    % Create a scatter plot
    figure('Position', [162, 233, 713, 592])
    s = scatterfit(MouseThalROICoords(:, XYZ), mpc_thal, 80, cmap(50,:), [], 0);
    
    % Calculate Pearson's correlation coefficient and related statistics
    [RHO, pval] = corr(MouseThalROICoords(:, XYZ), mpc_thal, 'Type', 'Pearson');
    [~, ~, CIL, CIU] = corrcoef(MouseThalROICoords(:, XYZ), mpc_thal);

    % Create figure caption
    FigCaption{XYZ} = [annotlabels{XYZ}, ', mPC', num2str(grad), ' score with ', ...
                       lower(SPATIAL_DIR{XYZ}), ' axis (CCFv3 ', SPATIAL_AXIS{XYZ}, ...
                       '-coordinate): Pearson''s r(', num2str(length(mpc_thal) - 2), ...
                       ') = ', num2str(RHO), ', p = ', num2str(pval), ...
                       ', CI = [', num2str(CIL(1, 2)), ', ', num2str(CIU(1, 2)), ...
                       '], two-tailed.'];
    
    % Display the caption
    disp(FigCaption{XYZ})

    % Set labels, titles, and annotations
    xlabel({SPATIAL_DIR{XYZ}, ['(CCFv3 {\it', SPATIAL_AXIS{XYZ}, '}-coordinate)']})
    ylabel(['Thalamic mouse PC', num2str(grad), ' score'])
    set(gca, 'FontSize', 24)
    xlimits = xlim;
    ylimits = ylim;
    text_x_coord = find_point_on_line(xlimits(1), xlimits(2), 0.05);
    text_y_coord = find_point_on_line(ylimits(1), ylimits(2), 1);

    % Format p-value for display
    if pval > .001
        p_val_format = ['{\itp} = ', num2str(round(pval, 3))];
    else
        p_rounded = num2str(pval, '%.2s');
        eloc = strfind(p_rounded, 'e');
        p_val_format = ['{\itp} = ', p_rounded(1:eloc - 1), '\times10^{', p_rounded(eloc + 1:end), '}'];
    end
    title(['{\itr} = ', num2str(round(RHO, 3)), ', ', p_val_format], 'FontWeight', 'Normal')
    a = annotation('textbox', [0, .896, .05, .13], 'String', annotlabels{XYZ}, ...
                   'FontSize', 32, 'EdgeColor', 'none');
    
    % Save figures as PNG and SVG
    print([SavePreFix, '_', annotlabels{XYZ}, '.png'], '-dpng', '-r300')
    print([SavePreFix, '_', annotlabels{XYZ}, '.svg'], '-dsvg')
end

end
