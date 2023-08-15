function [SigMapDescrips, SigMapDescripsFull] = PlotSigNeuromaps(neuromap_corrs, outdir, xlabel_var)
% PlotSigNeuromaps: A function to create scatter plots of significant
% neuromap correlations
% 
% Inputs:
%   neuromap_corrs - Struct containing neuromap correlation data
%   outdir (optional) - Directory to save the output figures (default: './figure_outputs/SigNeuromaps')
%   xlabel_var (optional) - Label for the x-axis variable (default: 'Cortical region PC1 loading')
%
% Outputs:
%   SigMapDescrips - Cell array containing short descriptions of significant neuromap correlations
%   SigMapDescripsFull - Cell array containing detailed descriptions of significant neuromap correlations
%

%   This function creates scatter plots of significant neuromap correlations
%   and saves them as PNG and SVG files. It also generates short and detailed descriptions
%   of the correlations, including correlation coefficients and p-values.

% Check if neuromap_corrs struct was provided
    if nargin < 1
        load('./data/processed/NeuroMapCorrs.mat', 'neuromap_corrs')
    end

    % Check if output directory was provided
    if nargin < 2
        outdir = './figure_outputs/SigNeuromaps';
    end

    % Check if xlabel variable was provided
    if nargin < 3
        xlabel_var = 'Cortical region PC1 loading';
    end
    
    % Create the output directory if it doesn't exist
    mkdir(outdir)

    % Find significant maps based on correlation data
    sig_maps = find(neuromap_corrs.corr_sig);

    Ndata = length(neuromap_corrs.data);
    
    NsigMaps = length(sig_maps);
    SigMapDescrips = cell(1,NsigMaps);
    SigMapDescripsFull = cell(1,NsigMaps);

    % Loop through significant maps and create scatter plots
    for i = 1:NsigMaps
        figure('Position', [291 391 560 420])
        SigInd = sig_maps(i);
        
        annotlabel = numberToLetter(i);
        
        if neuromap_corrs.p_perm(SigInd) == 0
            pval_format = '{\itp_{spin}} < .0001 ';   
        else
            pval_format = ['{\itp_{spin}} = ', num2str(round(neuromap_corrs.p_perm(SigInd), 4))];
        end
        s = scatterfit(neuromap_corrs.data, neuromap_corrs.neuromap_parc(:, SigInd), 40, neuromap_corrs.data, [], 0);
        title(['{\itr} = ', num2str(neuromap_corrs.corr(SigInd), 3), ', ', pval_format, ''], 'FontWeight', 'normal');
        xlabel(xlabel_var)
        ylabel(neuromap_corrs.name{SigInd})
        colormap(turbo)
        set(gca, 'FontSize', 16)
        s.MarkerEdgeColor = [0 0 0];
        a = annotation('textbox', [0 .905 .05 .13], 'String', annotlabel, 'FontSize', 32, 'EdgeColor', 'none');
        
        % Save the figures as PNG and SVG files
        print([outdir, '/Sigmap_', annotlabel, '.png'], '-dpng', '-r300')
        print([outdir, '/Sigmap_', annotlabel, '.svg'], '-dsvg')
        close all
        
        % Collect short and detailed descriptions of correlations
        SigMapDescrips{i} = neuromap_corrs.description{SigInd};
        
        [~, ~, CIL, CIU] = corrcoef(neuromap_corrs.data, neuromap_corrs.neuromap_parc(:, SigInd), 'Rows', 'complete');
        
        SigMapDescripsFull{i} = [annotlabel, ', ', neuromap_corrs.description{SigInd}, ' (Pearson''s r(', num2str(Ndata - 2), ...
            ') = ', num2str(neuromap_corrs.corr(SigInd), 3), ', pspin = ', num2str(round(neuromap_corrs.p_perm(SigInd), 4)), ...
            ', CI = [', num2str(CIL(1, 2), 3), ', ', num2str(CIU(1, 2), 3), '], two-tailed)'];  
    end
end
