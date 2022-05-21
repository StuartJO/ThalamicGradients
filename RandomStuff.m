

for i = 1:size(MZ_ID,1)
    IDs = MZ_ID(i,:);
    IDs_IND = find(ismember(IDs,TWINSUBs));
    for j = 1:2
        ID = MZ_ID(i,j);
        for k = 1:921
        if ~isnan(ID)
           ID_IND = find(TWINSUBs==ID);
           vec = Pall2{k}(ID_IND,:);
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
           vec = Pall2{k}(ID_IND,:);
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