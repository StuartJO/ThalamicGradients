function s = scatterfit(x,y,sz,c,pformat,plotCorr,plotCorrSig)

if nargin < 3
    sz = 36;
end

if nargin < 4
    c = [0 0.4470 0.7410];
end

if nargin < 5
    pformat = [];
end

if nargin < 6
    plotCorr = 1;
end

if nargin < 7
    plotCorrSig = 0;
end

s = scatter(x,y,sz,c,'filled');

[r,p] = corr(x,y,'rows','pairwise');

nan_inds = logical(isnan(x)+isnan(y));

% Get coefficients of a line fit through the data.
coefficients = polyfit(x(~nan_inds), y(~nan_inds),1);
% Create a new x axis with exactly 1000 points (or whatever you want).
xFit = linspace(nanmin(x), nanmax(x), 1000);
% Get the estimated yFit value for each of those 1000 new x locations.
yFit = polyval(coefficients , xFit);
% Plot everything.
hold on; % Set hold on so the next plot does not blow away the one we just drew.
plot(xFit, yFit, 'k-', 'LineWidth', 2); % Plot fitted line.

xlimits = xlim;
ylimits = ylim;
set(gca,'FontSize',18)
xlim(xlimits)
ylim(ylimits)

if isempty(pformat)

if p > .001
p_val_format = ['{\itp} = ',num2str(round(p,3))];
else
    if plotCorrSig == 0
p_val_format = '{\itp} < .001';    
    else
        p_rounded = num2str(p,'%.2s');        
        eloc = strfind(p_rounded,'e');    
        p_val_format = ['{\itp} = ',p_rounded(1:eloc-1),'\times10^{',p_rounded(eloc+1:end),'}'];            
    end
end

else
    p_val_format = pformat;
end

if r < 0
    d = .7;
else
    d = .1;
end
if plotCorr == 1
text_x_coord = find_point_on_line(xlimits(1),xlimits(2),d);
text_y_coord = find_point_on_line(ylimits(1),ylimits(2),.9);
text(text_x_coord,text_y_coord,{['{\itr} = ',num2str(round(r,3))],p_val_format},'FontSize',18)
elseif plotCorr == 2
    title(['{\itr} = ',num2str(round(r,3)),', ',p_val_format],'FontWeight','Normal');
end