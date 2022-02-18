
for i = 1:84
    data = ThalSeed(:,1:250,i);
    sub_data_norm{i} = zscore(data);
    [sub_coeff{i},sub_score{i},~,~,sub_explained{i}] = pca(sub_data_norm{i});
end

all_scores = [score{3} sub_score];

score_align = procrustes_alignment(sub_score,'reference',score{3});

C_align = zeros(85);

C_align(1,1) = 1;
for i = 1:84
    C_align(i+1,1) = corr(score_align{i}(:,1),score{3}(:,1));
    C_align(1,i+1) = C_align(i+1,1);
    for j = 1:84
    c = corr(score_align{i}(:,1),score_align{j}(:,1)); 
    C_align(i+1,j+1) = c; 
    C_align(j+1,i+1) = c;
    end
end


C_align = zeros(250,84);
for i = 1:84
    C_align(:,i) = corr(sub_score{i},score{3}(:,1));
end
[~,I] = max(abs(C_align));