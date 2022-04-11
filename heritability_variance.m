for i = 1:Ntwinsubs
    d = diag(cov(score_align{i}));
    score_align_explained(i) = d(1)/sum(d);
    d = diag(cov(score_align_250{i}));
    score_align_250_explained(i) = d(1)/sum(d);
    d = diag(cov(score_align_matlab{i}));
    score_align_matlab_explained(i) = d(1)/sum(d);
    d = diag(cov(score_align_matlab_obl{i}));
    score_align_matlab_obl_explained(i) = d(1)/sum(d);
    d = diag(cov(score_align_tract{i}));
    score_align_tract_explained(i) = d(1)/sum(d);
end

for i = 1:Ntwinsubs
for j = 1:3
pc_thal = zscore(sub_score{i}(:,j));
RHO(i,j) = corr(seed_mni_coords(:,1),pc_thal,'Type','Spearman');
end
[maxRHO(i),maxRHOIND(i)] = max(abs(RHO(i,:)));
end

load('twinCovariatesDWI.mat')
Output_MZ = nan(size(MZ_ID,1),2,7);
for i = 1:size(MZ_ID,1)
    for j = 1:2
        ID = MZ_ID(i,j);
        if ~isnan(ID)
           ID_IND = find(TWINSUBs==ID);
           Output_MZ(i,j,1) = score_align_explained(ID_IND); 
           Output_MZ(i,j,2) = score_align_250_explained(ID_IND); 
           Output_MZ(i,j,3) = score_align_matlab_explained(ID_IND); 
           Output_MZ(i,j,4) = score_align_matlab_obl_explained(ID_IND); 
           Output_MZ(i,j,5) = sub_explained{ID_IND}(1); 
           Output_MZ(i,j,6) = sub_explained{ID_IND}(2);  
           Output_MZ(i,j,7) = score_align_tract_explained(ID_IND); 
           Output_MZ(i,j,8) = maxRHO(ID_IND); 
           Output_MZ(i,j,9) = sub_explained{ID_IND}(maxRHOIND(ID_IND)); 
        end
    end
end

Output_DZ = nan(size(DZ_ID,1),2,7);
for i = 1:size(DZ_ID,1)
    for j = 1:2
        ID = DZ_ID(i,j);
        if ~isnan(ID)
           ID_IND = find(TWINSUBs==ID);
           Output_DZ(i,j,1) = score_align_explained(ID_IND); 
           Output_DZ(i,j,2) = score_align_250_explained(ID_IND); 
           Output_DZ(i,j,3) = score_align_matlab_explained(ID_IND); 
           Output_DZ(i,j,4) = score_align_matlab_obl_explained(ID_IND); 
           Output_DZ(i,j,5) = sub_explained{ID_IND}(1); 
           Output_DZ(i,j,6) = sub_explained{ID_IND}(2); 
           Output_DZ(i,j,7) = score_align_tract_explained(ID_IND);  
           Output_DZ(i,j,8) = maxRHO(ID_IND); 
           Output_DZ(i,j,9) = sub_explained{ID_IND}(maxRHOIND(ID_IND)); 
        end
    end
end

save('TwinAlignment_new.mat','Output_MZ','Output_DZ')

for i = 1:9
   subplot(3,3,i)
   scatter(Output_DZ(:,1,i),Output_DZ(:,2,i),10,'b','filled')
   hold on
   scatter(Output_MZ(:,1,i),Output_MZ(:,2,i),10,'r','filled')
   
   x = Output_DZ(:,1,i); y = Output_DZ(:,2,i);
   % Get coefficients of a line fit through the data.
coefficients = polyfit(x, y, 1);
% Create a new x axis with exactly 1000 points (or whatever you want).
xFit = linspace(min(x), max(x), 1000);
% Get the estimated yFit value for each of those 1000 new x locations.
yFit = polyval(coefficients , xFit);
% Plot everything.
hold on; % Set hold on so the next plot does not blow away the one we just drew.
plot(xFit, yFit, 'b', 'LineWidth', 2); % Plot fitted line.
grid on;

   x = Output_MZ(:,1,i); y = Output_MZ(:,2,i);
   % Get coefficients of a line fit through the data.
coefficients = polyfit(x, y, 1);
% Create a new x axis with exactly 1000 points (or whatever you want).
xFit = linspace(min(x), max(x), 1000);
% Get the estimated yFit value for each of those 1000 new x locations.
yFit = polyval(coefficients , xFit);
% Plot everything.
hold on; % Set hold on so the next plot does not blow away the one we just drew.
plot(xFit, yFit, 'r', 'LineWidth', 2); % Plot fitted line.
grid on;
end