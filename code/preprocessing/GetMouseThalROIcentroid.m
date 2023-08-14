load('./data/ancillary/MouseOhParc.mat')

ThalRegions = 88:122;

MouseThalOnly = MouseOhParc;
MouseThalOnly(~ismember(MouseOhParc,ThalRegions)) = NaN;
MouseThalOnly(MouseThalOnly==0) = NaN;
MouseThalOnly(1:228,:,:) = NaN;

MouseThalROI = changem(MouseThalOnly,1:35,ThalRegions);

MouseThalROICoords = zeros(35,3);
%
for i = 1:35
IND = find(MouseThalROI==i);
[mX,mY,mZ] = ind2sub(size(MouseThalROI),IND);
% mouse_medpos(i) = mean(mX);
d = squareform(pdist([mX,mY,mZ]));
[~,I] = min(mean(d));
MouseThalROICoords(i,:) = [mX(I) mY(I) mZ(I)];
disp(num2str(i))
end

save('./data/ancillary/MouseThalROICoords.mat','MouseThalROICoords')