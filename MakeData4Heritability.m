
Output_MZ = nan(size(MZ_ID,1),4,2);
for i = 1:size(MZ_ID,1)
    for j = 1:4
        ID = MZ_ID(i,j);
        if ~isnan(ID)
           ID_IND = find(TWINSUBs==ID)+1;
           Output_MZ(i,j,1) = C_align(ID_IND,1); 
           Output_MZ(i,j,2) = C_align(ID_IND,1); 
        end
    end
end

Output_DZ = nan(size(DZ_ID,1),4,2);
for i = 1:size(DZ_ID,1)
    for j = 1:4
        ID = DZ_ID(i,j);
        if ~isnan(ID)
           ID_IND = find(TWINSUBs==ID)+1;
           Output_DZ(i,j,1) = C_align(ID_IND,1); 
           Output_DZ(i,j,2) = C_align(ID_IND,1); 
        end
    end
end

save('TwinAlignment.mat','Output_MZ','Output_DZ')