
for i = 1:494
data = TwinThalSeed(:,1:250,i);
norm_data = BF_NormalizeMatrix(data,'scaledSigmoid');
%norm_data = data;
norm_data(isnan(norm_data))=0;
for j = 1:921
twin_seed_con{j}(i,:) = norm_data(j,:);
end
end

%Pall{i} = squareform(pdist(twin_vec));
for j = 1:921
Pall{j} = corr(twin_seed_con{j}');
Pall_1{j} = (corr(twin_seed_con{j}')+1)./2;
Pall_2{j} = squareform(pdist(twin_seed_con{j}));

end

thrs = [0 .25 .5 .75];
for t = 1:4
for j = 1:921
d = twin_seed_con{j};
for k = 1:494
    for l = 1:494
        A = double(d(k,:)>thrs(t));
        B = double(d(l,:)>thrs(t));
        DC(k,l) = 2*(sum(A.*B))/sum(A + B);
    end
end
P_dc{t}{j} = DC;
end
end

for t = 1:4
for j =1:921
K = GenSimilarity(Kind,Kind);
P = P_dc{t}{j}(inGenSim,inGenSim);
X = covar(inGenSim,:);

[h2_dc{t}(j)] = h2_mat(P, K, X,0);

end
end

GenSimilarity = csvread('K.csv');

GenSimilaritySub = GenSimilarity(1,:);
GenSimilarity(1,:) = [];

count = 1;
inGenSim = true(494,1);

age = [DZ_age(:); MZ_age(:)];
sex = [DZ_sex(:); MZ_sex(:)];
covID = [DZ_ID(:); MZ_ID(:)];

for i = 1:494

ind = find(GenSimilaritySub==TWINSUBs(i));
if isempty(ind)
    inGenSim(i) = false;
else
    Kind(count) = ind;
    count = count + 1;
end
covar_ind = find(covID==TWINSUBs(i));
covar(i,:) = [age(covar_ind(1)) sex(covar_ind(1))];
end

for j =1:921
K = GenSimilarity(Kind,Kind);
P = Pall{j}(inGenSim,inGenSim);
X = covar(inGenSim,:);

[h2(j), p_perm(j),] = h2_mat(P, K, X,0);

P = Pall_1{j}(inGenSim,inGenSim);

[h2_1(j), p_perm_1(j),] = h2_mat(P, K, X,0);

P = Pall_2{j}(inGenSim,inGenSim);

[h2_2(j), p_perm_1(j),] = h2_mat(P, K, X,0);
end


load('MNI_Seed_voxelData.mat','seeds_vox','seed_coords')

seed_voxel_coords = seeds_vox(logical(seed_ind),:);
PlotThalGradient3(h2,seed_voxel_coords,turbo(256),'h2',2.1)

h2_parc = h2;
h2_parc(h2<.5) = 1;
h2_parc(h2>=.5) = 2;

norm_data = BF_NormalizeMatrix(ThalSeedAvg,'scaledSigmoid');

parc1_norm_data = norm_data(:,1:250).*repmat(double(h2_parc'==1),1,250);

parc2_norm_data = norm_data(:,1:250).*repmat(double(h2_parc'==2),1,250);

parc1_norm_data(parc1_norm_data==0) = NaN;
parc2_norm_data(parc2_norm_data==0) = NaN;
parc1_mean = nanmean(parc1_norm_data)*-1;
parc2_mean = nanmean(parc2_norm_data);

parc_weights = [parc1_mean; parc2_mean];

[parc_h2,parc_assign_ind] = max([abs(parc1_mean); parc2_mean]);

parc_h2(parc_assign_ind==1) = parc_h2(parc_assign_ind==1)*-1;

load('fsaverage_surface_data.mat')
surface.vertices = lh_inflated_verts;
surface.faces = lh_faces;
cmap = turbo(256);
    
figure('Position',[461   462   560   325])
p = plotSurfaceROIBoundary(surface,lh_rand500,parc_h2,'midpoint',cmap,1,2,[-1 1]);
camlight(80,-10);
camlight(-80,-10);
view([-90 0])
colorbar

axis off
axis image


load('twinCovariatesDWI.mat')
for g = 1:2
Thal_Output_MZ = nan(size(MZ_ID,1),2,921);
for i = 1:size(MZ_ID,1)
    for j = 1:2
        ID = MZ_ID(i,j);
        if ~isnan(ID)
           ID_IND = find(TWINSUBs==ID);
           Thal_Output_MZ(i,j,:) = zscore(sub_score{ID_IND}(:,g)); 
        end
    end
end

Thal_Output_DZ = nan(size(DZ_ID,1),2,921);
for i = 1:size(DZ_ID,1)
    for j = 1:2
        ID = DZ_ID(i,j);
        if ~isnan(ID)
           ID_IND = find(TWINSUBs==ID);
           Thal_Output_DZ(i,j,:) = zscore(sub_score{ID_IND}(:,g)); 
        end
    end
end


for i = 1:921
   
   x = Thal_Output_DZ(:,1,i); y = Thal_Output_DZ(:,2,i);

DZ_pcorr = partialcorr(x,y,[DZ_sex(:,1) DZ_age(:,1)]);
DZ_corr = corr(x,y);

   x = Thal_Output_MZ(:,1,i); y = Thal_Output_MZ(:,2,i);

MZ_pcorr = partialcorr(x,y,[MZ_sex(:,1) MZ_age(:,1)]);
MZ_corr = corr(x,y);

Thal_h(i,g) = 2*(MZ_corr-DZ_corr);
Thal_h_partial(i,g) = 2*(MZ_pcorr-DZ_pcorr);

end

end


Pall{j} = corr(twin_seed_con{j}');
Pall_1{j} = (corr(twin_seed_con{j}')+1)./2;
Pall_2{j} = squareform(pdist(twin_seed_con{j}));


for i = 1:size(MZ_ID,1)
    IDs = MZ_ID(i,:);
    IDs_IND = find(ismember(IDs,TWINSUBs));
    for j = 1:2
        ID = MZ_ID(i,j);
        for k = 1:921
        if ~isnan(ID)
           ID_IND = find(TWINSUBs==ID);
           vec = Pall{k}(ID_IND,:);
           vec(IDs_IND) = NaN;
           Pall_Output_MZ(i,j,k) = nanmean(vec); 
        end
        end
    end
end

for i = 1:size(DZ_ID,1)
    IDs = DZ_ID(i,:);
    IDs_IND = find(ismember(IDs,TWINSUBs));
    for j = 1:2
        ID = DZ_ID(i,j);
        for k = 1:921
        if ~isnan(ID)
           ID_IND = find(TWINSUBs==ID);
           vec = Pall{k}(ID_IND,:);
           vec(IDs_IND) = NaN;
           Pall_Output_DZ(i,j,k) = nanmean(vec); 
        end
        end
    end
end

for i = 1:921
     
   x = Pall_Output_DZ(:,1,i); y = Pall_Output_DZ(:,2,i);

DZ_pcorr = partialcorr(x,y,[DZ_sex(:,1) DZ_age(:,1)]);
DZ_corr = corr(x,y);

   x = Pall_Output_MZ(:,1,i); y = Pall_Output_MZ(:,2,i);

MZ_pcorr = partialcorr(x,y,[MZ_sex(:,1) MZ_age(:,1)]);
MZ_corr = corr(x,y);

Pall_h(i) = 2*(MZ_corr-DZ_corr);
Pall_h_partial(i) = 2*(MZ_pcorr-DZ_pcorr);

end





for i = 1:size(MZ_ID,1)
    IDs = MZ_ID(i,:);
    IDs_IND = find(ismember(IDs,TWINSUBs));
    for j = 1:2
        ID = MZ_ID(i,j);
        for k = 1:921
        if ~isnan(ID)
           ID_IND = find(TWINSUBs==ID);
           vec = Pall_2{k}(ID_IND,:);
           vec(IDs_IND) = NaN;
           Pall2_Output_MZ(i,j,k) = nanmean(vec); 
        end
        end
    end
end

for i = 1:size(DZ_ID,1)
    IDs = DZ_ID(i,:);
    IDs_IND = find(ismember(IDs,TWINSUBs));
    for j = 1:2
        ID = DZ_ID(i,j);
        for k = 1:921
        if ~isnan(ID)
           ID_IND = find(TWINSUBs==ID);
           vec = Pall_2{k}(ID_IND,:);
           vec(IDs_IND) = NaN;
           Pall2_Output_DZ(i,j,k) = nanmean(vec); 
        end
        end
    end
end

for i = 1:921
     
   x = Pall2_Output_DZ(:,1,i); y = Pall2_Output_DZ(:,2,i);

DZ_pcorr = partialcorr(x,y,[DZ_sex(:,1) DZ_age(:,1)]);
DZ_corr = corr(x,y);

   x = Pall2_Output_MZ(:,1,i); y = Pall2_Output_MZ(:,2,i);

MZ_pcorr = partialcorr(x,y,[MZ_sex(:,1) MZ_age(:,1)]);
MZ_corr = corr(x,y);

Pall2_h(i) = 2*(MZ_corr-DZ_corr);
Pall2_h_partial(i) = 2*(MZ_pcorr-DZ_pcorr);

end


seed_voxel_coords = seeds_vox(logical(seed_ind),:);

for i = 1:2
PlotThalGradient3(Thal_h_partial(:,i),seed_voxel_coords,turbo(256),'h2',2.1)
print(['./NewFigures/Tract_PC',num2str(i),'_h2.png'],'-dpng')
figure
histogram(Thal_h_partial(:,i))
ylabel('h2')
print(['./NewFigures/Tract_PC',num2str(i),'_hist_h2.png'],'-dpng')
end



for t = 1:4
PlotThalGradient3(h2_dc{t},seed_voxel_coords,turbo(256),'h2',2.1)
print(['./NewFigures/DC_thr',num2str(t),'_h2.png'],'-dpng')
figure
histogram(h2_dc{t})
ylabel('h2')
print(['./NewFigures/DC_thr',num2str(t),'_hist_h2.png'],'-dpng')
end





seed_voxel_coords = seeds_vox(logical(seed_ind),:);
PlotThalGradient3(h2,seed_voxel_coords,turbo(256),'h2',2.1)
print(['./NewFigures/Corr_h2.png'],'-dpng')
figure
histogram(h2)
ylabel('h2')
print(['./NewFigures/Corr_hist_h2.png'],'-dpng')
h2_parc = h2;
h2_parc(h2<.5) = 1;
h2_parc(h2>=.5) = 2;

norm_data = BF_NormalizeMatrix(ThalSeedAvg,'scaledSigmoid');

parc1_norm_data = norm_data(:,1:250).*repmat(double(h2_parc'==1),1,250);

parc2_norm_data = norm_data(:,1:250).*repmat(double(h2_parc'==2),1,250);

parc1_norm_data(parc1_norm_data==0) = NaN;
parc2_norm_data(parc2_norm_data==0) = NaN;
parc1_mean = nanmean(parc1_norm_data)*-1;
parc2_mean = nanmean(parc2_norm_data);

parc_weights = [parc1_mean; parc2_mean];

[parc_h2,parc_assign_ind] = max([abs(parc1_mean); parc2_mean]);

parc_h2(parc_assign_ind==1) = parc_h2(parc_assign_ind==1)*-1;

load('fsaverage_surface_data.mat')
surface.vertices = lh_inflated_verts;
surface.faces = lh_faces;
cmap = turbo(256);
    
figure('Position',[461   462   560   325])
p = plotSurfaceROIBoundary(surface,lh_rand500,parc_h2,'midpoint',cmap,1,2,[-1 1]);
camlight(80,-10);
camlight(-80,-10);
view([-90 0])
colorbar

axis off
axis image

print(['./NewFigures/Corr_h2_cort.png'],'-dpng')