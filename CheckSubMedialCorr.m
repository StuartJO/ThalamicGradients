
for i = 1:Ntwinsubs
for j = 1:3
pc_thal = zscore(sub_score{i}(:,j));
RHO(i,j) = corr(seed_mni_coords(:,1),pc_thal,'Type','Spearman');
end
[maxRHO(i),maxRHOIND(i)] = max(abs(RHO(i,:)));
end

