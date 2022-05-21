
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

score_align = procrustes_alignment(sub_score,'reference',score{3});



for i = 1:100
    data = TwinThalSeed(:,1:250,i);
    for j = 1
    A = Bipart_rewire_cols(data,10);
    norm_data = BF_NormalizeMatrix(A,'scaledSigmoid');
    norm_data(isnan(norm_data)) = 0;
    Rewired{i,j} = norm_data;
    [Rewired_coeff{i,j},Rewired_score{i,j},~,~,Rewired_explained{i,j}] = pca(norm_data);
    end
    disp(['Finished ',num2str(i),'/',num2str(Ntwinsubs)])
end

score_align_rewired = procrustes_alignment(Rewired_score(1:100,1),'reference',score{3});

for i = 1:100
score_align_rewired_cov = cov(score_align_rewired{i});
score_align_rewired_explained(i,:) = (diag(score_align_rewired_cov)./sum(diag(score_align_rewired_cov)))*100;
end

score_align = procrustes_alignment(sub_score(1:100),'reference',score{3});

for i = 1:100
score_align_rewired_matlab{i} = rotatefactors(Rewired_score{i,1},'Method','procrustes','Target',score{3}(:,1:250),'Type','orthogonal');
end

for i = 1:100
score_align_rewired_cov = cov(score_align_rewired_matlab{i});
score_align_rewired_matlab_explained(i,:) = (diag(score_align_rewired_cov)./sum(diag(score_align_rewired_cov)))*100;
end



for i = 1:100
score_align_cov = cov(score_align{i});
score_align_explained(i,:) = (diag(score_align_cov)./sum(diag(score_align_cov)))*100;
end

% 
% for i = 1:size(MZ_ID,1)
%     for j = 1:2
%         ID = MZ_ID(i,j);
%            ID_IND = find(TWINSUBs==ID);
%     data = TwinThalSeed(:,1:250,ID_IND);
%     for iter = 1
%     A = Bipart_rewire_cols(data,10);
%     norm_data = BF_NormalizeMatrix(data,'scaledSigmoid');
%     norm_data(isnan(norm_data)) = 0;
%     Rewired{i,iter} = norm_data;
%     [Rewired_coeff{i,iter},Rewired_score{i,iter},~,~,Rewired_explained{i,iter}] = pca(norm_data);
%     end
%     disp(['Finished ',num2str(i),'/',num2str(Ntwinsubs)])
%         end
% end