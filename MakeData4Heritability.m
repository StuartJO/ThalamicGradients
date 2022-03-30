load('twinCovariatesDWI.mat')
Output_MZ = nan(size(MZ_ID,1),4,3);
for i = 1:size(MZ_ID,1)
    for j = 1:4
        ID = MZ_ID(i,j);
        if ~isnan(ID)
           ID_IND = find(TWINSUBs==ID)+1;
           Output_MZ(i,j,1) = C_align(ID_IND,1); 
           Output_MZ(i,j,2) = abs(C_orig(1,ID_IND-1)); 
           Output_MZ(i,j,3) = abs(C_orig(2,ID_IND-1));
        end
    end
end

Output_DZ = nan(size(DZ_ID,1),4,3);
for i = 1:size(DZ_ID,1)
    for j = 1:4
        ID = DZ_ID(i,j);
        if ~isnan(ID)
           ID_IND = find(TWINSUBs==ID)+1;
           Output_DZ(i,j,1) = C_align(ID_IND,1); 
           Output_DZ(i,j,2) = abs(C_orig(1,ID_IND-1)); 
           Output_DZ(i,j,3) = abs(C_orig(2,ID_IND-1));
        end
    end
end

save('TwinAlignment.mat','Output_MZ','Output_DZ')