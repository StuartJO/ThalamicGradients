#!/bin/env bash

module load mrtrix
module load fsl

DIR="/projects/kg98/stuarto/SeedReg/dHCP/NewSub"
WORKDIR="/projects/kg98/stuarto/SeedReg/dHCP/NewSub"
MM=1.5;
RERUN=0;

SUB="CC00179XX15"
SES="58800"

if [ $RERUN = "1" ]; then

fslmaths ${DIR}/anat/sub-${SUB}_ses-${SES}_desc-drawem87_dseg.nii.gz -thr 43 -uthr 43 -bin ${DIR}/Left_Thal.nii

flirt -in ${DIR}/Left_Thal.nii.gz -ref ${DIR}/Left_Thal.nii.gz -applyisoxfm ${MM} -nosearch -out ${DIR}/Left_Thal_${MM}mm.nii.gz
fslmaths ${DIR}/Left_Thal_${MM}mm.nii.gz -thr .5 -bin ${DIR}/Left_Thal_${MM}mm_masked.nii.gz

cp ${DIR}/anat/sub-${SUB}_ses-${SES}_desc-drawem9_dseg.nii.gz ${DIR}/drawem9_dseg.nii.gz

for i in {1..9}; do
fslmaths ${DIR}/drawem9_dseg.nii.gz -thr ${i} -uthr ${i} -bin ${DIR}/seg_${i}.nii.gz
done
fslmaths ${DIR}/seg_1.nii.gz -add ./seg_5.nii.gz -bin ${DIR}/csf.nii.gz
fslmaths ${DIR}/seg_2.nii.gz -add ./seg_9.nii.gz -add ${DIR}/seg_6.nii.gz -bin ${DIR}/gm.nii.gz
fslmaths ${DIR}/seg_3.nii.gz -add ./seg_8.nii.gz -bin ${DIR}/wm.nii.gz
fslmaths ${DIR}/seg_7.nii.gz -bin ${DIR}/gm_sub.nii.gz
fslmaths ${DIR}/drawem9_dseg.nii.gz -mul 0 ${DIR}/path.nii.gz

mrcat ${DIR}/gm.nii.gz ${DIR}/gm_sub.nii.gz ${DIR}/wm.nii.gz ${DIR}/csf.nii.gz ${DIR}/path.nii.gz ${DIR}/ACT.nii

tckgen -algorithm Seedtest -seed_grid_per_voxel ${DIR}/Left_Thal_${MM}mm_masked.nii.gz 1 -output_seeds ${DIR}/seeds_${MM}mm.txt ${DIR}/drawem9_dseg.nii.gz ${DIR}/seeds_${MM}mm_tracts.tck -force -nthreads 0

nlines=$(wc -l ${DIR}/seeds_${MM}mm.txt | awk '{ print $1 }')
nseeds="$(($nlines - 2))"
tck_ind=$(seq 1 1 ${nseeds})
echo $tck_ind > ${DIR}/seeds_${MM}mm_ind.txt

fi

SEEDS="${WORKDIR}/seeds_1.5mm.txt"

mkdir ${WORKDIR}/tracts/

nlines=$(wc -l ${SEEDS} | awk '{ print $1 }')
nseeds="$(($nlines - 2))"

START=1
END="$nseeds"
echo $END
for (( c=$START; c<=$END; c++ ))
do
seedid="$(($c+1))"
seeddata=$(sed -n "${seedid}p" ${SEEDS})
X="$(cut -d',' -f3 <<<"${seeddata}")"
Y="$(cut -d',' -f4 <<<"${seeddata}")"
Z="$(cut -d',' -f5 <<<"${seeddata}")"

#tckgen ${DIR}/FOD.mif ./seedtck_10000/thal_seed_${c}.tck -seed_sphere $X,$Y,$Z,0 -select 10000 -seed_unidirectional -act ${DIR}/ACT.nii -backtrack -include ./Cortex.nii
#echo $X,$Y,$Z
tckgen ${WORKDIR}/ss_wmfod_eddy.mif ${WORKDIR}/tracts/thal_seed_${c}.tck -seed_sphere $X,$Y,$Z,0 -select 5000 -seed_unidirectional -act ${WORKDIR}/ACT.nii -backtrack

#tckgen ${WORKDIR}/ss_wmfod_eddy.mif ${WORKDIR}/tracts/thal_seed_${c}.tck -seed_sphere $X,$Y,$Z,0 -select 5000 -seed_unidirectional -include gm.nii.gz 

#tckgen ${WORKDIR}/ss_wmfod_eddy.mif ${WORKDIR}/tracts/thal_seed_${c}.tck -seed_sphere -7.344,16.14,20.81,0 -select 5000 -seed_unidirectional -act ${WORKDIR}/ACT.nii -backtrack
#tckgen ${WORKDIR}/ss_wmfod_eddy.mif ${WORKDIR}/test.tck -seed_image ${WORKDIR}/wm.nii.gz -exclude csf.nii.gz -select 1000000 -include gm.nii.gz 
#tckgen ${WORKDIR}/ss_wmfod_eddy.mif ${WORKDIR}/test.tck -seed_image ${WORKDIR}/gm_sub.nii.gz -select 1000000 -act ${WORKDIR}/ACT.nii -backtrack

#tck2connectome ${WORKDIR}/tracts/thal_seed_${c}.tck $PARC ${WORKDIR}/tracts/thal_seed_${c}_cp -vector -assignment_radial_search 5 -force

done

#msmapplywarp /projects/kg98/stuarto/SeedReg/dHCP/week-40_hemi-left_space-dhcpSym_dens-32k_sphere.surf.gii ${DIR}/week-40-${SUB}_wm -anat ${DIR}/xfm/sub-${SUB}_ses-${SES}_hemi-left_from-native_to-dhcpSym40_dens-32k_mode-sphere.surf.gii ${DIR}/anat/sub-${SUB}_ses-${SES}_hemi-left_wm.surf.gii 

# Loop over tck files and add them to a string

RUN=0

if [ $RUN = "1" ]; then 

THAL="./tracts/thal_seed_1.tck"
for (( i=2; i<=$END; i++ ))
do 
THAL="${THAL} ./tracts/thal_seed_${i}.tck"
done

# Combine all the tck files. Doing *tck will do the same but streamlines are not ordered. nthreads needs to be 0 for streamlines to be ordered as well

tckedit ${THAL} thal_seed.tck -nthreads 0

SeedStart=1
for (( c=$START; c<=$END; c++ ))
do
Nstream=$(tckstats -output count thal_seed_${c}.tck -quiet)
SeedEnd=$(($SeedStart + $Nstream - 1))

if [ $c -eq 1 ]; then
echo "$SeedStart $SeedEnd $Nstream" > SeedStreamIndex
else
echo "$SeedStart $SeedEnd $Nstream" >> SeedStreamIndex
fi

SeedStart=$(($SeedEnd + 1))

done

fi
