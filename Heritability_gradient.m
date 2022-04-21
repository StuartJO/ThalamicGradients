
load('twinCovariatesDWI.mat')
Thal_Output_MZ = nan(size(MZ_ID,1),2,921);
for i = 1:size(MZ_ID,1)
    for j = 1:2
        ID = MZ_ID(i,j);
        if ~isnan(ID)
           ID_IND = find(TWINSUBs==ID);
           Thal_Output_MZ(i,j,:) = zscore(score_align{ID_IND}(:,1)); 
        end
    end
end

Thal_Output_DZ = nan(size(DZ_ID,1),2,921);
for i = 1:size(DZ_ID,1)
    for j = 1:2
        ID = DZ_ID(i,j);
        if ~isnan(ID)
           ID_IND = find(TWINSUBs==ID);
           Thal_Output_DZ(i,j,:) = zscore(score_align{ID_IND}(:,1)); 
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

Thal_h(i) = 2*(MZ_corr-DZ_corr);
Thal_h_partial(i) = 2*(MZ_pcorr-DZ_pcorr);

end

PlotThalGradient3((score{3}(:,1)),seed_voxel_coords,turbo(256),['Heritability'],2.1)
figure('Position',[162   233   713   592])
seed_mni_coords = seed_coords(logical(seed_ind))';
SPATIAL_DIR = {'Medial-Lateral','Anterior-posterior','Dorsal-ventral'};
      pc_thal = zscore(score{PCtype}(:,grad));
    s = scatter(seed_mni_coords(:,XYZ),pc_thal,'filled','MarkerFaceAlpha',.25,'MarkerFaceColor',cmap(50,:));
    [RHO,pval] = corr(seed_mni_coords(:,XYZ),pc_thal,'Type','Spearman');
    %title([SPATIAL_DIR{i},' ',data_type_name{j},' \rho = ',num2str(RHO)])
    xlabel({'Medial-lateral position','(MNI x-axis coordinate)'})
    ylabel(['Thalamic seed PC',num2str(grad),' score'])
    set(gca,'FontSize',24,'XDir','reverse')
    text(-0.5,1.5,['{\itr_{s}} = ',num2str(round(RHO,3))],'FontSize',24)
    pval_string = num2str(pval,3);
    
    pval_m = pval_string(1:find(pval_string=='e')-1);
    if length(pval_m) < 4
        pval_m = [pval_m,'0'];
    end
    pval_n = pval_string(find(pval_string=='e')+1:end);
    text(-0.5,1,['{\itp} = ',pval_m,'\times10^{',pval_n,'}'],'FontSize',24)
