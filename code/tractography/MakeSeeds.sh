#!/bin/env bash

module load mrtrix
module load fsl
module load ants/2.2.0 

DIR="/projects/kg98/stuarto/SeedReg"

SUBJECT_LIST="/projects/kg98/stuarto/SeedReg/Unrelated_HCP_subs_100.txt"

RERUN=0

MM=1.75

if [ $RERUN = "1" ]; then

fslmaths ${DIR}/Tian_Subcortex_S1_3T_1mm.nii.gz -thr 11 -uthr 12 -bin ${DIR}/Left_Thal.nii
cp /usr/local/fsl/6.0.4/fsl/data/standard/MNI152_T1_1mm_brain.nii.gz ${DIR}/brain_1mm.nii.gz
cp /usr/local/fsl/6.0.4/fsl/data/standard/MNI152_T1_2mm_brain.nii.gz ${DIR}/brain_2mm.nii.gz

# Make ${MM} seeds

flirt -in ${DIR}/Left_Thal.nii.gz -ref ${DIR}/Left_Thal.nii.gz -applyisoxfm ${MM} -nosearch -out ${DIR}/Left_Thal_${MM}mm.nii.gz
fslmaths ${DIR}/Left_Thal_${MM}mm.nii.gz -thr .5 -bin ${DIR}/Left_Thal_${MM}mm_masked.nii.gz

tckgen -algorithm Seedtest -seed_grid_per_voxel ${DIR}/Left_Thal_${MM}mm_masked.nii.gz 1 -output_seeds ${DIR}/seeds_${MM}mm.txt ${DIR}/brain_1mm.nii.gz ${DIR}/seeds_${MM}mm_tracts.tck -force -nthreads 0

nlines=$(wc -l ${DIR}/seeds_${MM}mm.txt | awk '{ print $1 }')
nseeds="$(($nlines - 2))"
tck_ind=$(seq 1 1 ${nseeds})
echo $tck_ind > ${DIR}/seeds_${MM}mm_ind.txt

fi

nsubs=$(wc -l ${SUBJECT_LIST} | awk '{ print $1 }')

# FYI, "fsl" is in the filename because at one stage I was exploring using ANTS to do the registration and so named the files appropriately to keep track. I stuck with FSL because ANTS appeared to distort things too much

for ID in $(seq 1 $nsubs); do
SUB=$(sed -n "${ID}p" ${SUBJECT_LIST})
WORKDIR="${DIR}/$SUB"
mkdir ${WORKDIR}
warpconvert /projects/hcp1200_processed/2021/Preprocessed/${SUB}/MNINonLinear/xfms/acpc_dc2standard.nii.gz displacement2deformation ${WORKDIR}/warp_MNIto${SUB}.nii.gz -quiet -force
fslmaths /projects/hcp1200_processed/2021/Preprocessed/${SUB}/T1w/parc/aparc+first_acpc.nii -thr 35 -uthr 35 -bin ${WORKDIR}/${SUB}_ThalMask.nii.gz

tcktransform ${DIR}/seeds_${MM}mm_tracts.tck ${WORKDIR}/warp_MNIto${SUB}.nii.gz ${WORKDIR}/${SUB}_seeds_${MM}mm_fsl.tck -force -quiet -nthreads 0 
tckedit -include ${WORKDIR}/${SUB}_ThalMask.nii.gz -tck_weights_in ${DIR}/seeds_${MM}mm_ind.txt -tck_weights_out ${WORKDIR}/masked_seed_ind_${MM}mm_fsl_${SUB}.txt ${WORKDIR}/${SUB}_seeds_${MM}mm_fsl.tck ${WORKDIR}/seeds_masked_${MM}mm_fsl_${SUB}.tck -quiet -force -nthreads 0
tr -s ' '  '\n'< ${WORKDIR}/masked_seed_ind_${MM}mm_fsl_${SUB}.txt > ${WORKDIR}/masked_seed_ind_${MM}mm_fsl_rows_${SUB}.txt
done

tckedit -include ${DIR}/GeneMask.nii.gz -tck_weights_in ${DIR}/seeds_${MM}mm_ind.txt -tck_weights_out ${DIR}/gene_masked_seeds_${MM}mm_ind.txt ${DIR}/seeds_${MM}mm_tracts.tck ${DIR}/seeds_gene_masked_${MM}mm_tracts.tck

#tckedit -tck_weights_in ${DIR}/seeds_${MM}mm_fsl_seeds2use_valid_bin.txt ${WORKDIR}/${SUB}_seeds_${MM}mm_fsl.tck ${WORKDIR}/seeds_masked_${MM}mm_fsl_${SUB}.tck -quiet -force -nthreads 0
