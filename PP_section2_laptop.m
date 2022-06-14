load('Sub76_ThalData.mat')

figure('Position',[0 0 1920 1080])

EntrezIDs_all = dlmread('AHBAEntrez.txt');

GeneIDs_all_table = readtable('AHBAgeneSymbol.txt','ReadVariableNames',false);

GeneIDs_all = GeneIDs_all_table.Var1;

[~,GeneIDs_all_Gryglewski_ind] = ismember(GenesKept,EntrezIDs_all);

GeneIDs_all_valid = GeneIDs_all(GeneIDs_all_Gryglewski_ind);

MouseGeneTable = readtable('MouseThalGenes.xlsx');

MouseGenes = MouseGeneTable.GeneSymbol;

MouseGenesUpper = upper(MouseGenes);

idx = find(ismember(GeneIDs_all_valid, MouseGenesUpper));

MouseGenesValid = GeneIDs_all_valid(idx);

% data_type{1} = zscore(ThalSeedAvg(:,1:250));
% data_type{2} = zscore(SeedGene_kept);

data_type{1} = BF_NormalizeMatrix(ThalSeedAvg(:,1:250),'scaledSigmoid');
data_type{2} = BF_NormalizeMatrix(SeedGene_kept,'scaledSigmoid');

data_type{3} = [data_type{1} data_type{2}];
data_type{4} = data_type{2}(:,idx);
data_type{5} = [data_type{1} data_type{2}(:,idx)];

load('fsaverage_surface_data.mat')
load('lh.avgMapping_allSub_RF_ANTs_MNI152_orig_to_fsaverage.mat');
lh_mni_vox = MNI_mm2vox(ras','mm');
load('rh.avgMapping_allSub_RF_ANTs_MNI152_orig_to_fsaverage.mat');
rh_mni_vox = MNI_mm2vox(ras','mm');
MNIsurface.vertices = lh_mni_vox;
MNIsurface.faces = lh_faces;

MNIsurfaceR.vertices = rh_mni_vox;
MNIsurfaceR.faces = rh_faces;
load('MNI_Seed_voxelData.mat')
brain = niidata;

pMNIsurface = patch(MNIsurface);
set(pMNIsurface,'EdgeColor','none','FaceColor',[.5 .5 .5],'Clipping','off');
pMNIsurface.FaceLighting = 'gouraud';
material dull

view([180 0])

axis off
axis tight
axis equal
axis vis3d

hold on
pMNIsurfaceR = patch(MNIsurfaceR);
set(pMNIsurfaceR,'EdgeColor','none','FaceColor',[.5 .5 .5],'Clipping','off');
pMNIsurfaceR.FaceLighting = 'gouraud';
material dull

clight = camlight('headlight');

xlimits = xlim;
ylimits = ylim;
zlimits = zlim;

end_xlimits = [133.3088  143.8052];
end_ylimits = [133.1712  162.0478];
end_zlimits = [81.3852  105.1673];

slice_xlimits = [98.0211  108.5175];
slice_ylimits = [93.4533  127.2028];
slice_zlimits = [116.0609  143.8561];


t= 60;
views = [linspace(180,0,t)' linspace(0,90,t)'];
xlimits_range = [linspace(xlimits(1),slice_xlimits(1),t)' linspace(xlimits(2),slice_xlimits(2),t)'];
ylimits_range = [linspace(ylimits(1),slice_ylimits(1),t)' linspace(ylimits(2),slice_ylimits(2),t)'];
zlimits_range = [linspace(zlimits(1),slice_zlimits(1),t)' linspace(zlimits(2),slice_zlimits(2),t)'];


TR = stlread('LeftThalMeshSmoothManual.stl');
Thalsurf.faces = TR.ConnectivityList;
Thalsurf.vertices = MNI_mm2vox(TR.Points,'mm');

% Thalsurf.faces = TR.faces;
% Thalsurf.vertices = MNI_mm2vox(TR.vertices,'mm');


pThalsurf = patch(Thalsurf);
set(pThalsurf,'EdgeColor','none','FaceColor',[1 .647 0],'Clipping','off');
pThalsurf.FaceLighting = 'gouraud';
material dull

seed_voxel_coords = seeds_vox(logical(seed_ind),:);
z_loc = 82;

load('MNI_1mm_mask2.mat','brain_mask')

brain_rescaled = rescale(brain,1,256);
brain_rescaled(~brain_mask)=NaN;
brain_slice = brain_rescaled(:,:,z_loc)';
colormap(gray(256))
slicesurf = surf(repmat(z_loc, size(brain_slice)), brain_slice);
slicesurf.EdgeAlpha=0;
slicesurf.FaceAlpha=0;
slicesurf.FaceLighting='none';
slicesurf.Clipping='off';

xlim(xlimits)
ylim(ylimits)
zlim(zlimits)


Cortalphavals = [1 1 1 1 1 1 1 1 1 1 linspace(1,0,50)];
CortalphavalsR = [1 1 1 1 1 1 1 1 1 1 linspace(1,0,50)];
%Thalalphavals = [0 0 0 0 0 0 0 0 0 0 ones(1,50)];
slicealphavals = [zeros(1,50) linspace(0,1,10)];
for i = 1:60
    pThalsurf.FaceAlpha = 0;
    %pThalsurf.FaceAlpha = Thalalphavals(i);
    pMNIsurface.FaceAlpha = Cortalphavals(i);
    pMNIsurfaceR.FaceAlpha = Cortalphavals(i);
    view(views(i,:))
    xlim(xlimits_range(i,:));
    ylim(ylimits_range(i,:));
    zlim(zlimits_range(i,:));
    pause(.1)    
    delete(clight)
    slicesurf.FaceAlpha=slicealphavals(i);
    clight = camlight('headlight');
    pause(.1)
    print(['./GIF/PP2_S1_',num2str(i),'.png'],'-dpng')
end



%%
voxel_seeds = seeds_vox(logical(seed_ind),:);
orange_cmap=repmat([1 .647 0],256,1);

z_locs = [82 82 83:90 89:-1:69 70:82];
for j = 1:2
    if exist('sSeeds','var')
    delete(sSeeds)
    end
dist2slice = pdist2(seed_voxel_coords(:,3),z_locs(j));
minDist = min(dist2slice);
seeds2plot_ind = dist2slice==minDist;

[x,y,z]=sphere(64);
sphere_mesh = surf2patch(x,y,z,'triangles');

seed_voxel_coords = seeds_vox(logical(seed_ind),:);

seeds2plot = seed_voxel_coords(seeds2plot_ind,:);

[newslicevals,newCmap] = ProjThalData2Slice(ones(921,1),z_locs(j),voxel_seeds,orange_cmap,2);
slicesurf.CData=newslicevals;
colormap(newCmap)
caxis([1 512])

if j > 1 
    clear sSeeds
for i = 1:length(seeds2plot)
    sphere_mesh_roi = sphere_mesh;
    sphere_mesh_roi.vertices = (sphere_mesh.vertices.*.25) + [seeds2plot(i,1:2) 83];
    sSeeds(i) = patch(sphere_mesh_roi,'EdgeColor','none','FaceColor',[0 0 0]./255,'Clipping','off','FaceLighting','gouraud');
    material dull     
end
end
pause(.1)
print(['./GIF/PP2_S2_',num2str(j),'.png'],'-dpng')

end

%%

gene_data = SeedGene_kept;
gene_data = data_type{2};

genes2use = [1:30 round(linspace(31,2233,91))];
voxel_seeds = seeds_vox(logical(seed_ind),:);
thalCmap=viridis(256);
for i = 1:120
thalData = zscore(gene_data(:,genes2use(i)));   
[newslicevals,newCmap] = ProjThalData2Slice(thalData,z_loc,voxel_seeds,thalCmap,2);
colormap(newCmap)
slicesurf.CData=newslicevals;
caxis([1 512])
pause(.1)
print(['./GIF/PP2_S3_',num2str(i),'.png'],'-dpng')
end

slicesurf.CData=brain_rescaled(:,:,z_loc)';
colormap(gray(256))
caxis([1 256])
print(['./GIF/PP2_S4.png'],'-dpng')
delete(clight)
camlight(80,-10);
camlight(-80,-10);

xlimits_side = [22.8127  156.9987];
ylimits_side = [24.5933  195.4591];
zlimits_side = [27.7807  149.4851];

t = 60;
views = [linspace(0,90,t)' linspace(90,0,t)'];
xlimits_range = [linspace(slice_xlimits(1),xlimits_side(1),t)' linspace(slice_xlimits(2),xlimits_side(2),t)'];
ylimits_range = [linspace(slice_ylimits(1),ylimits_side(1),t)' linspace(slice_ylimits(2),ylimits_side(2),t)'];
zlimits_range = [linspace(slice_zlimits(1),zlimits_side(1),t)' linspace(slice_zlimits(2),zlimits_side(2),t)'];

Slicealphavals_rev = [0 0 0 0 0 0 0 0 0 0 linspace(0,1,50)];

alphavals = [0 0 0 0 0 0 0 0 0 0 linspace(0,.25,50)];

tracks = read_mrtrix_tracks('192035inMNI_thal_seed_620_cortical.tck');
clear tracts_MNI
for i = 1:length(tracks.data)
tracts_MNI{i} = MNI_mm2vox(tracks.data{i},'mm');
end


seed2plot = seed_voxel_coords(620,:);
sphere_mesh_roi = sphere_mesh;
sphere_mesh_roi.vertices = (sphere_mesh.vertices.*.25) + [seed2plot(1,1:2)+.5 83];
RSeeds = patch(sphere_mesh_roi,'EdgeColor','none','FaceColor',[255 0 0]./255,'Clipping','off','FaceLighting','none');
material dull   
    print(['./GIF/PP2_S5.png'],'-dpng')
    
for j = 1:t
    if exist('streamline_handle','var')
    delete(streamline_handle)
    end
    for i = 1:length(tracts_MNI)
        s = tracts_MNI{i};
        streamline_length = size(s,1);

        Sr = round(find_point_on_line(1,streamline_length,j/60));
        
        S = s(1:Sr,:)+[.5 .5 0];
        
        if size(S,1) > 1

        [~,SDIR] = max(sum(abs(diff(S))));    

        if SDIR == 1
            SCOL = [1 0 0 .2];
            col = 'r';
        elseif SDIR == 2
            SCOL = [0 1 0 .2];
            col = 'g';
        else
            SCOL = 	[0 0 1 .2];
            col = 'b';
        end

        StepRate = 1;

        streamline_handle(i) = patch([S(:,1)' NaN],[S(:,2)' NaN],[S(:,3)' NaN],0); 
        joint = ([S(2:end,:); NaN NaN NaN]-S);
        joint(end,:) = joint(end-1,:);
        temp_joint = joint;
        joint(:,1) = temp_joint(:,1);
        joint(:,2) = temp_joint(:,2);
        cdata = [abs(joint./StepRate); NaN NaN NaN];
        cdata = reshape(cdata,length(cdata),1,3);
        set(streamline_handle(i),'CData', cdata, 'EdgeColor','interp','FaceColor','interp','Clipping','off') 

        hold on

        end
    end

    for k = 1:length(sSeeds)
    sSeeds(k).FaceAlpha=1-Slicealphavals_rev(j);
    end
    
slicesurf.FaceAlpha=1-Slicealphavals_rev(j);
pMNIsurface.FaceAlpha = alphavals(j);
pMNIsurfaceR.FaceAlpha = alphavals(j);
    view(views(j,:))
    xlim(xlimits_range(j,:));
    ylim(ylimits_range(j,:));
    zlim(zlimits_range(j,:));
    pause(.1)
    print(['./GIF/PP2_S6_',num2str(j),'.png'],'-dpng')
end

delete(streamline_handle)

conn_data = (ThalSeedAvg(620,1:250));

conn_data = norm(620,:);

FaceVertexCData = makeFaceVertexCData(MNIsurface.vertices,MNIsurface.faces,lh_rand500,conn_data,parula(256));
BOUNDARY = findROIboundaries(MNIsurface.vertices,MNIsurface.faces,lh_rand500,'midpoint');
pMNIsurface.FaceVertexCData = FaceVertexCData;
pMNIsurface.FaceColor='interp';
pMNIsurface.FaceAlpha=1;
for i = 1:length(BOUNDARY)
   boundary_plot.boundary(i) = plot3(BOUNDARY{i}(:,1), BOUNDARY{i}(:,2), BOUNDARY{i}(:,3), 'Color', 'k', 'LineWidth',2,'Clipping','off');
end
print(['./GIF/PP2_S7.png'],'-dpng')

rowspace = [2:30 round(linspace(31,921,91))];

for i = 1:length(rowspace)
    conn_data = norm(rowspace(i),:);
FaceVertexCData = makeFaceVertexCData(MNIsurface.vertices,MNIsurface.faces,lh_rand500,conn_data,parula(256));
pMNIsurface.FaceVertexCData = FaceVertexCData;    
print(['./GIF/PP2_S8_',num2str(i),'.png'],'-dpng')
end

figure

imagesc(norm(seed_order,tract_order))
xticks([])
yticks([])
colormap(parula(1000))

order_tract_data = tract_data(seed_order,tract_order);
order_tract_data_orig=order_tract_data;


order_tract_data(2:end,:) = NaN;
tract_plot = pcolor(order_tract_data);
tract_plot.EdgeAlpha=0;
set(gca,'YDir','reverse')
% xlim([1 10])
% ylim([-3 6])

xlinspace = [10 10 linspace(10,251,118)];

ylinspace = [xlinspace(1:30)-4 linspace(xlinspace(31)-4,921,90)];

ylinspace2 = linspace(-3,1,120);

rowspace = [1:30 round(linspace(31,921,90))];
% yaxisscale = (1:90).^10;
% ylinspace = rescale(yaxisscale./max(yaxisscale),6,921);
axis square
%axis off
xticks([])
yticks([])

% set(gca,'Color','w')
% 
% f=getframe(gca);
% imwrite(f.cdata,['./GIF/PP2_S9.png'],'png');
set(gca,'XColor','none')
set(gca,'YColor','none')

%exportgraphics(gca,['./GIF/PP2_S9.png'])

f=getframe(gca);
imwrite(f.cdata,['./GIF/PP2_S9.png'],'png');

for i = 1:120
    tract_plot.CData(1:rowspace(i),:) = order_tract_data_orig(1:rowspace(i),:);
%     xlim([1 xlinspace(i)])
%     ylim([ylinspace2(i) ylinspace(i)])
    pause(.1)
    f=getframe(gca);
imwrite(f.cdata,['./GIF/PP2_S10_',num2str(i),'.png'],'png');
    %exportgraphics(gcf,['./GIF/PP2_S10_',num2str(i),'.png'])
end
% 
% ax1 = gca;
% 
% P = get(ax1,'pos');
% axes_mat = axes('Position',P);
% set(axes_mat,'color','none')

genes2use = [1:30 round(linspace(31,2233,91))];

figure
gene_data = data_type{2};
imagesc(gene_data(seed_order,gene_order))
xticks([])
yticks([])
colormap(viridis(1000))

order_gene_data = gene_data(seed_order,gene_order);
order_gene_data_orig=order_gene_data;


order_gene_data(:,2:end) = NaN;
tract_plot = pcolor(order_gene_data);
tract_plot.EdgeAlpha=0;
set(gca,'YDir','reverse')
xlim([1 10])
ylim([1 10])

xlinspace = linspace(10,2233,120);

ylinspace = linspace(10,921,120);

ylinspace2 = linspace(-3,1,120);

rowspace = [round(linspace(1,2233,120))];
% yaxisscale = (1:90).^10;
% ylinspace = rescale(yaxisscale./max(yaxisscale),6,921);
axis square
%axis off
xticks([])
yticks([])

% set(gca,'Color','w')
% 
% f=getframe(gca);
% imwrite(f.cdata,['./GIF/PP2_S9.png'],'png');
set(gca,'XColor','none')
set(gca,'YColor','none')

%exportgraphics(gca,['./GIF/PP2_S9.png'])

f=getframe(gca);
imwrite(f.cdata,['./GIF/GENE_MAT.png'],'png');

for i = 1:120
    tract_plot.CData(:,1:rowspace(i)) = order_gene_data_orig(:,1:rowspace(i));
    xlim([1 xlinspace(i)])
    ylim([1 ylinspace(i)])
    pause(.1)
    f=getframe(gca);
imwrite(f.cdata,['./GIF/GENE_MAT_',num2str(i),'.png'],'png');
    %exportgraphics(gcf,['./GIF/PP2_S10_',num2str(i),'.png'])
end




f=getframe(gca);
imwrite(f.cdata,['./GIF/GENE_MAT.png'],'png');

for i = 1:120
    tract_plot.CData(:,1:rowspace(i)) = order_gene_data_orig(:,1:rowspace(i));
    xlim([1 2233])
    ylim([1 921])
    pause(.1)
    f=getframe(gca);
imwrite(f.cdata,['./GIF/GENE_MAT_',num2str(i),'.png'],'png');
    %exportgraphics(gcf,['./GIF/PP2_S10_',num2str(i),'.png'])
end