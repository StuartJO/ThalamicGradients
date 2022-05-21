function [slicedata,newCmap] = ProjThalData2Slice(thalData,z_loc,voxel_seeds,thalCmap,thr)

load('MNI_Seed_voxelData.mat','thalmask','niidata')

inmask = find(thalmask==1);
[X,Y,Z] = ind2sub(size(thalmask),inmask);
ThalMaskDists = pdist2([X Y Z],voxel_seeds);
[ThalMask_nearest_seedmm,ThalMask_nearest_seed] = min(ThalMaskDists');
ThalVoxClust = changem(ThalMask_nearest_seed,thalData,1:length(voxel_seeds));
thalmaskclust = thalmask;
thalmaskclust(thalmaskclust==0) = NaN;
thalmaskclust(inmask) = ThalVoxClust;
thalmaskclust(inmask(ThalMask_nearest_seedmm>thr)) = NaN;

CmapSize = size(thalCmap,1);

newCmap = [gray(CmapSize); thalCmap];

thalmaskclust_scaled = rescale(thalmaskclust,CmapSize+1,CmapSize*2);

thalslice = thalmaskclust_scaled(:,:,z_loc)';

brain_scaled = rescale(niidata,1,CmapSize);

slicedata = brain_scaled(:,:,z_loc)';

slicedata(find(~isnan(thalslice))) = thalslice(find(~isnan(thalslice)));
