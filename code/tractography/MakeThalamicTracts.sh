#!/bin/env bash
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH -t 1-0:0:0
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --export=ALL
#SBATCH --mem-per-cpu=12G
#SBATCH --account=kg98
#SBATCH --output=/projects/kg98/stuarto/ThalamicGradients/SLURM_OUTPUT/slurm-%j.out

# Note you mnay need to make the SLURM_OUTPUT directory
module load mrtrix
module load fsl

## Below is how you run in in bash

#SUBJECT_LIST="/projects/kg98/stuarto/ThalamicGradients/VALIDSEED_UnrelatedSubs.txt"
#nsubs=$(wc -l ${SUBJECT_LIST} | awk '{ print $1 }')
#for ID in $(seq 1 $nsubs); do SUB=$(sed -n "${ID}p" ${SUBJECT_LIST}); sbatch ./MakeThalamicTracts.sh $SUB; done

SUB=$1

WORKDIR="/projects/kg98/stuarto/ThalamicGradients/data/tractography/SUBJECTS/${SUB}"
SEEDS="${WORKDIR}/${SUB}_921seeds_1.75mm.txt"

if [ ! -d "${WORKDIR}/tracts_921seeds/" ]; then

#rm -r ${WORKDIR}/tracts/
mkdir ${WORKDIR}/tracts_921seeds

nseeds=$(wc -l ${SEEDS} | awk '{ print $1 }')
#nseeds="$(($nlines - 2))"

START=1
END="$nseeds"

PARC="/projects/hcp1200_processed/2021/Preprocessed/${SUB}/T1w/parc/random500_acpc.nii"
PARENTDIR="/projects/hcp1200_processed/2021/Processed/${SUB}"
 
for (( c=$START; c<=$END; c++ ))
do
seedid="$(($c))"
seeddata=$(sed -n "${seedid}p" ${SEEDS})
#X="$(cut -d',' -f3 <<<"${seeddata}")"
#Y="$(cut -d',' -f4 <<<"${seeddata}")"
#Z="$(cut -d',' -f5 <<<"${seeddata}")"

#tckgen ${DIR}/FOD.mif ./seedtck_10000/thal_seed_${c}.tck -seed_sphere $X,$Y,$Z,0 -select 10000 -seed_unidirectional -act ${DIR}/ACT.nii -backtrack -include ./Cortex.nii
tckgen ${PARENTDIR}/FOD.mif ${WORKDIR}/tracts_921seeds/thal_seed_${c}.tck -seed_sphere ${seeddata},0 -select 5000 -seed_unidirectional -act ${PARENTDIR}/ACT.nii -backtrack

tck2connectome ${WORKDIR}/tracts_921seeds/thal_seed_${c}.tck $PARC ${WORKDIR}/tracts_921seeds/thal_seed_${c}_cp -vector -assignment_radial_search 5 -force

done

fi
