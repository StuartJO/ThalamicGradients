function PlotTractogramPCs()

load('./data/processed/main_decomp.mat','pcs_thal','pcs_cort')

numColors = 256;
ramp = linspace(-100,100, numColors);
cform = makecform('lab2srgb');
a = repmat(ramp, [numColors 1]);           % -a on left
b = repmat(flipud(ramp'), [1 numColors]);  % -b on bottom
L = 50 * ones(numColors, numColors);  % A single L value.
Lab = cat(3, L, a, b); % A 2D image.
colormap2D = applycform(Lab, cform);

PC_THAL = pcs_thal(:,1);
PC_CORT = pcs_cort(:,1);

score_cmap_rescaled = rescale(zscore(PC_THAL),1,256);

coeff_cmap_rescaled = rescale(zscore(PC_CORT),1,256);

figure

imagesc(colormap2D);
set(gca,'YDir','normal')
xlabel('Thalamic seed PC1 score')
ylabel('Cortical ROI PC1 coefficient')

xticks_rescaled = rescale([-2:1 min(zscore(PC_THAL)) max(zscore(PC_THAL))],1,256);
yticks_rescaled = rescale([-2:1 min(zscore(PC_CORT)) max(zscore(PC_CORT))],1,256);

xticks(xticks_rescaled(1:end-2))
yticks(yticks_rescaled(1:end-2))

set(gca,'XTickLabel',-2:1)
set(gca,'YTickLabel',-2:1)
set(gca,'FontSize',20)
axis square

print('./figure_outputs/TractogramCmap.svg','-dsvg','-r300')
print('./figure_outputs/TractogramCmap.png','-dpng','-r300')
% 

figure('Position',[0 0 1920 1080])

camlight(80,-10);
camlight(-80,-10);
view([90 0])

axis off
axis tight
axis equal
axis vis3d

f = waitbar(0,'Starting');
MNIExampleTracks = cell(921,3);
totalTractLengths = 0;

for j = 1:921
    tracks = read_mrtrix_tracks (['./data/ExampleTracks/MNI/seed_',num2str(j),'.tck']);
    ROI = dlmread(['./data/ExampleTracks/conn/seed_',num2str(j),'_assignments']);
    tracks_data = tracks.data;
    Ntracks = length(tracks_data);
    Tracklen = zeros(length(tracks_data),1);
    rmv = false(Ntracks,1);
    
    for i = 1:Ntracks   
        Tracklen(i) = size(tracks_data{i},1);
        if ~(ROI(i) <= 250 && ROI(i) ~= 0)
            rmv(i) = true;        
        end
    end
    
    tracks_data(rmv) = [];
    ROI(rmv) = [];
    Tracklen(rmv) = [];
    MNIExampleTracks{j,1} = tracks_data;
    MNIExampleTracks{j,2} = ROI;
    MNIExampleTracks{j,3} = Tracklen;
    
    totalTractLengths = totalTractLengths + sum(Tracklen) + length(Tracklen);
    
    waitbar(j/921,f,['Finished seed ',num2str(j),' of ',num2str(921)]);
end
close(f)

STREAMLINES = zeros(totalTractLengths,3);
STREAMLINES_COLOR = zeros(totalTractLengths,3);

f = waitbar(0,'Starting');
IND = 0;
for j = 1:921
    tracks = MNIExampleTracks{j,1};
    ROI = MNIExampleTracks{j,2};
    
    for i = 1:length(tracks)
            tract = tracks{i};
            TractLength = length(tract)+IND;
            color2use = squeeze(colormap2D(round(score_cmap_rescaled(j)),round(coeff_cmap_rescaled(ROI(i))),:))';

            STREAMLINES(IND+1:TractLength+1,:) = [[tract(:,1)' NaN]',[tract(:,2)' NaN]',[tract(:,3)' NaN]'];
            STREAMLINES_COLOR(IND+1:TractLength+1,:) = repmat(color2use,size(tract,1)+1,1);

            IND = TractLength + 1;
    end
    waitbar(j/921,f,['Finished seed ',num2str(j),' of ',num2str(921)]);
end
close(f)

cdata_all = reshape(STREAMLINES_COLOR,length(STREAMLINES_COLOR),1,3);    
streamlines = patch(STREAMLINES(:,1),STREAMLINES(:,2),STREAMLINES(:,3),0);
set(streamlines,'CData', cdata_all, 'EdgeColor','interp','FaceColor','interp','Clipping','off') 

exportgraphics(gcf,'./figure_outputs/TractogramLateral.png','Resolution',300)

view([90 90])

exportgraphics(gcf,'./figure_outputs/TractogramSuperior.png','Resolution',300)