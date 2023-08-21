#!/bin/env bash

module load mrtrix
module load fsl
module load ants/2.2.0 
/projects/kg98/stuarto/ThalamicGradients/data/
TRACTDIR="tractography"

PREPROCESSEDDIR="/projects/kg98/stuarto/ThalamicGradients/data/preprocessed"

#mkdir ${TRACTDIR}

SUBJECT_LIST="/projects/kg98/stuarto/ThalamicGradients/UnrelatedSubs.txt"

RERUN=0

MM=1.75

if [ $RERUN = "1" ]; then

fslmaths ${PREPROCESSEDDIR}/Tian_Subcortex_S1_3T_1mm.nii.gz -thr 11 -uthr 12 -bin ${PREPROCESSEDDIR}/Left_Thal.nii

# Get the T1w's from FSL and copy them just for later ease of use
cp /usr/local/fsl/6.0.4/fsl/data/standard/MNI152_T1_1mm_brain.nii.gz ${TRACTDIR}/brain_1mm.nii.gz
cp /usr/local/fsl/6.0.4/fsl/data/standard/MNI152_T1_2mm_brain.nii.gz ${TRACTDIR}/brain_2mm.nii.gz

# Make ${MM} seeds

# Yeah I could/should have done nearest neighbour interpolation with flirt but I forgot to
flirt -in ${PREPROCESSEDDIR}/Left_Thal.nii.gz -ref ${PREPROCESSEDDIR}/Left_Thal.nii.gz -applyisoxfm ${MM} -nosearch -out ${PREPROCESSEDDIR}/Left_Thal_${MM}mm.nii.gz
fslmaths ${PREPROCESSEDDIR}/Left_Thal_${MM}mm.nii.gz -thr .5 -bin ${PREPROCESSEDDIR}/Left_Thal_${MM}mm_masked.nii.gz

tckgen -algorithm Seedtest -seed_grid_per_voxel ${PREPROCESSEDDIR}/Left_Thal_${MM}mm_masked.nii.gz 1 -output_seeds ${PREPROCESSEDDIR}/seeds_${MM}mm.txt ${PREPROCESSEDDIR}/brain_1mm.nii.gz ${TRACTDIR}/seeds_${MM}mm_tracts.tck -force -nthreads 0

nlines=$(wc -l ${PREPROCESSEDDIR}/seeds_${MM}mm.txt | awk '{ print $1 }')
nseeds="$(($nlines - 2))"
tck_ind=$(seq 1 1 ${nseeds})
echo $tck_ind > ${PREPROCESSEDDIR}/seeds_${MM}mm_ind.txt

fi

nsubs=$(wc -l ${SUBJECT_LIST} | awk '{ print $1 }')

# FYI, "fsl" is in the filename because at one stage I was exploring using ANTS to do the registration and so named the files appropriately to keep track. I stuck with FSL because ANTS appeared to distort things too much

for ID in $(seq 1 $nsubs); do
SUB=$(sed -n "${ID}p" ${SUBJECT_LIST})
WORKDIR="${TRACTDIR}/$SUB"
mkdir ${WORKDIR}
warpconvert /projects/hcp1200_processed/2021/Preprocessed/${SUB}/MNINonLinear/xfms/acpc_dc2standard.nii.gz displacement2deformation ${WORKDIR}/warp_MNIto${SUB}.nii.gz -quiet -force
# This is each subjects segemnted brain where the cortex was segmented with Freesurfer but the subcortex with FSL. We extract just the thalamus to use as a mask for QC
fslmaths /projects/hcp1200_processed/2021/Preprocessed/${SUB}/T1w/parc/aparc+first_acpc.nii -thr 35 -uthr 35 -bin ${WORKDIR}/${SUB}_ThalMask.nii.gz

tcktransform ${TRACTDIR}/seeds_${MM}mm_tracts.tck ${WORKDIR}/warp_MNIto${SUB}.nii.gz ${WORKDIR}/${SUB}_seeds_${MM}mm_fsl.tck -force -quiet -nthreads 0 
tckedit -include ${WORKDIR}/${SUB}_ThalMask.nii.gz -tck_weights_in ${TRACTDIR}/seeds_${MM}mm_ind.txt -tck_weights_out ${WORKDIR}/masked_seed_ind_${MM}mm_fsl_${SUB}.txt ${WORKDIR}/${SUB}_seeds_${MM}mm_fsl.tck ${WORKDIR}/seeds_masked_${MM}mm_fsl_${SUB}.tck -quiet -force -nthreads 0
tr -s ' '  '\n'< ${WORKDIR}/masked_seed_ind_${MM}mm_fsl_${SUB}.txt > ${WORKDIR}/masked_seed_ind_${MM}mm_fsl_rows_${SUB}.txt
done

tckedit -include ${PREPROCESSEDDIR}/GeneMask.nii.gz -tck_weights_in ${PREPROCESSEDDIR}/seeds_${MM}mm_ind.txt -tck_weights_out ${PREPROCESSEDDIR}/gene_masked_seeds_${MM}mm_ind.txt ${TRACTDIR}/seeds_${MM}mm_tracts.tck ${TRACTDIR}/seeds_gene_masked_${MM}mm_tracts.tck

#tckedit -tck_weights_in ${TRACTDIR}/seeds_${MM}mm_fsl_seeds2use_valid_bin.txt ${WORKDIR}/${SUB}_seeds_${MM}mm_fsl.tck ${WORKDIR}/seeds_masked_${MM}mm_fsl_${SUB}.tck -quiet -force -nthreads 0
