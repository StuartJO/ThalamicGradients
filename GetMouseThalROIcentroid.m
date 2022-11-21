MouseThalOnly = MouseAtlasNew;
MouseThalOnly(~ismember(MouseAtlasNew,mouse_P56_id(ThalRegions))) = NaN;
MouseThalOnly(MouseThalOnly==0) = NaN;
MouseThalOnly(1:228,:,:) = NaN;

ThalRegions = 88:122;
MouseThalROI = changem(MouseThalOnly,1:35,mouse_P56_id(ThalRegions));

for i = 1:35
IND = find(MouseThalROI==i);
[mX,mY,mZ] = ind2sub(size(MouseThalROI),IND);
% mouse_medpos(i) = mean(mX);
d = squareform(pdist([mX,mY,mZ]));
[~,I] = min(mean(d));
MouseThalROICoord(i,:) = [mX(I) mX(I) mX(I)];
i
end

save('MouseThalROICoord.mat','MouseThalROICoord')