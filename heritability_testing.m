
data = TwinThalSeed(:,1:250,i);
norm_data = BF_NormalizeMatrix(data,'scaledSigmoid');

for i = 1:size(MZ_ID,1)
    for j = 1:2
        ID = MZ_ID(i,j);
        if ~isnan(ID)
           ID_IND = find(TWINSUBs==ID)+1;
%            Output_MZ(i,j,1) = C_align(ID_IND,1); 
%            Output_MZ(i,j,2) = abs(C_orig(1,ID_IND-1)); 
%            Output_MZ(i,j,3) = abs(C_orig(2,ID_IND-1));
%           mz_data(i,j) = new_explained(ID_IND-1);
           mz_data(i,j) = sub_explained{ID_IND-1}(1);
        end
    end
end

for i = 1:size(DZ_ID,1)
    for j = 1:2
        ID = DZ_ID(i,j);
        if ~isnan(ID)
           ID_IND = find(TWINSUBs==ID)+1;
%            Output_DZ(i,j,1) = C_align(ID_IND,1); 
%            Output_DZ(i,j,2) = abs(C_orig(1,ID_IND-1)); 
%            Output_DZ(i,j,3) = abs(C_orig(2,ID_IND-1));
%           dz_data(i,j) = new_explained(ID_IND-1);
           dz_data(i,j) = sub_explained{ID_IND-1}(1);
        end
    end
end

histogram(diff(mz_data,1,2))

hold on

histogram(diff(dz_data,1,2))