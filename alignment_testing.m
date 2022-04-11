
%addpath(genpath('C:\Users\Stuart\Desktop\ThalTractSeed_data\BrainSpace-master'))
addpath(genpath('F:\Documents\ThalTractSeed\BrainSpace-master'))
load('Seed_voxelData.mat')
load('MNI_Seed_voxelData.mat')
load('Sub76_ThalData.mat')
%load('ThalTractGeneData_926seeds_84subs.mat')
%load('Sub_ThalSeed.mat')
load('TwinSub_ThalData.mat')

Ntwinsubs = size(TwinThalSeed,3);

for i = 1:Ntwinsubs
    data = TwinThalSeed(:,1:250,i);
    norm_data = BF_NormalizeMatrix(data,'scaledSigmoid');
    norm_data(isnan(norm_data)) = 0;
    sub_data_norm{i} = norm_data;
    [sub_coeff{i},sub_score{i},~,~,sub_explained{i}] = pca(sub_data_norm{i});
end

data_type{1} = BF_NormalizeMatrix(ThalSeedAvg(:,1:250),'scaledSigmoid');
data_type{2} = BF_NormalizeMatrix(SeedGene_kept,'scaledSigmoid');
data_type{3} = [data_type{1} data_type{2}];

for j = 1:3
[coeff{j},score{j},~,~,explained{j}] = pca(data_type{j});
end

all_scores = [score{3} sub_score];

score_align_tract = procrustes_alignment(sub_score,'reference',score{1});

score_align = procrustes_alignment(sub_score,'reference',score{3});

score_align_250 = procrustes_alignment(sub_score,'reference',score{3}(:,1:250));

for i = 1:Ntwinsubs
score_align_matlab{i} = rotatefactors(sub_score{i},'Method','procrustes','Target',score{3}(:,1:250),'Type','orthogonal');
end

for i = 1:Ntwinsubs
score_align_matlab_obl{i} = rotatefactors(sub_score{i},'Method','procrustes','Target',score{3}(:,1:250));
end

for i = 1:Ntwinsubs
    d = diag(cov(score_align{i}));
    new_explained(i) = d(1)/sum(d);
    new_explained2(i) = d(3)/sum(d);
end


C_align = zeros(Ntwinsubs+1);

C_align(1,1) = 1;
for i = 1:Ntwinsubs
    C_align(i+1,1) = corr(score_align{i}(:,1),score{3}(:,1),'Type','Spearman');
    C_align(1,i+1) = C_align(i+1,1);
    for j = 1:Ntwinsubs
    c = corr(score_align{i}(:,1),score_align{j}(:,1),'Type','Spearman'); 
    C_align(i+1,j+1) = c; 
    C_align(j+1,i+1) = c;
    end
end


C_orig = zeros(250,Ntwinsubs);
for i = 1:Ntwinsubs
    C_orig(:,i) = corr(sub_score{i},score{3}(:,1),'Type','Spearman');
end
[~,I] = max(abs(C_orig));

figure('Position',[348 436 1462 664]);
for i = 1:10
   subplot(2,5,i)
   scatter(score{3}(:,1),sub_score{i}(:,1))
   xlabel('Group PC1')
   ylabel('Subject PC1')
   set(gca,'FontSize',16)
end
print(['./Figures/UnalignedPC1scatter.png'],'-dpng')

figure('Position',[348 436 1462 664]);
for i = 1:10
   subplot(2,5,i)
   scatter(score{3}(:,1),sub_score{i}(:,2))
   xlabel('Group PC1')
   ylabel('Subject PC2')
   set(gca,'FontSize',16)
end
print(['./Figures/UnalignedPC2scatter.png'],'-dpng')

figure('Position',[348 436 1462 664]);
for i = 1:10
   subplot(2,5,i)
   scatter(score{3}(:,1),score_align{i}(:,1))
   xlabel('Group PC1')
   ylabel('Subject aligned PC1')
   set(gca,'FontSize',16)
end
print(['./Figures/AlignedPC1scatter.png'],'-dpng')

load('MNI_Seed_voxelData.mat','seeds_vox')

seed_voxel_coords = seeds_vox(logical(seed_ind),:);

PlotThalGradient3(score{3}(:,3),seed_voxel_coords,turbo(100),'Group PC1',2.1)
print(['./Figures/GroupPC1.png'],'-dpng')

for i = 1:10
PlotThalGradient3(sub_score{i}(:,1),seed_voxel_coords,turbo(100),['Subject ',num2str(i),' unaligned PC1'],2.1)
print(['./Figures/Subject',num2str(i),'_unalignedPC1.png'],'-dpng')
end

for i = 1:10
PlotThalGradient3(score_align{i}(:,1),seed_voxel_coords,turbo(100),['Subject ',num2str(i),' aligned PC1'],2.1)
print(['./Figures/Subject',num2str(i),'_alignedPC1.png'],'-dpng')
end

figure('Position',[348 436 1176 420]);
subplot(1,3,1)
histogram(abs(C_orig(1,:)))
xlabel({'Subject unaligned PC1','correlation with Group PC1'})
set(gca,'FontSize',16)
subplot(1,3,2)
histogram(abs(C_orig(2,:)))
xlabel({'Subject unaligned PC2','correlation with Group PC1'})
set(gca,'FontSize',16)
subplot(1,3,3)
histogram(abs(C_align(1,2:end)))
xlabel({'Subject aligned PC1','correlation with Group PC1'})
set(gca,'FontSize',16)
print(['./Figures/CorrelationHistograms.png'],'-dpng')