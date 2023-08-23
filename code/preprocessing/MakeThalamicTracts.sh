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
#SBATCH --output=./SLURM_OUTPUT/slurm-%j.out

# Note you may need to make the SLURM_OUTPUT directory
module load mrtrix
module load fsl

## Below is how you run in in bash

#SUBJECT_LIST="/projects/kg98/stuarto/ThalamicGradients/VALIDSEED_UnrelatedSubs.txt"
#nsubs=$(wc -l ${SUBJECT_LIST} | awk '{ print $1 }')
#for ID in $(seq 1 $nsubs); do SUB=$(sed -n "${ID}p" ${SUBJECT_LIST}); sbatch ./MakeThalamicTracts.sh $SUB; done

SUB=$1

WORKDIR="./data/tractography/SUBJECTS/${SUB}"
SEEDS="${WORKDIR}/${SUB}_921seeds_1.75mm.txt"

if [ ! -d "${WORKDIR}/tracts_921seeds/" ]; then

#rm -r ${WORKDIR}/tracts/
mkdir ${WORKDIR}/tracts_921seeds
mkdir ${WORKDIR}/tracts_921seeds/random500
mkdir ${WORKDIR}/tracts_921seeds/Schaefer400_17net

nseeds=$(wc -l ${SEEDS} | awk '{ print $1 }')
#nseeds="$(($nlines - 2))"

START=1
END="$nseeds"

PARCDIR="/projects/hcp1200_processed/2021/Preprocessed/${SUB}/T1w/parc"
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

	for PARC in random500 Schaefer400_17net; do

		PARCPATH="${PARCDIR}/${PARC}_acpc.nii"
		tck2connectome ${WORKDIR}/tracts_921seeds/thal_seed_${c}.tck $PARCPATH ${WORKDIR}/tracts_921seeds/${PARC}/thal_seed_${c}_cp -vector -assignment_radial_search 5 -force

	done

done

fi

# Ok so obviously this is a complete waste of space, regenerating tracts for seeds which we have already generated them for. However, in defence, I am lazy. This came about because I just didn't want to figure out the indexing/mapping. In hindsight it is an easy fix, but ¯\_(ツ)_/¯

WORKDIR="./data/tractography/SUBJECTS/${SUB}"
SEEDS="${WORKDIR}/${SUB}_ALLGENEseeds_1.75mm.txt"

mkdir ${WORKDIR}/tracts_ALLGENEseeds
mkdir ${WORKDIR}/tracts_ALLGENEseeds/random500

if [ ! -d "${WORKDIR}ALLGENEseeds/" ]; then

mkdir ${WORKDIR}

nseeds=$(wc -l ${SEEDS} | awk '{ print $1 }')

START=1
END="$nseeds"

PARCDIR="/projects/hcp1200_processed/2021/Preprocessed/${SUB}/T1w/parc"
PARENTDIR="/projects/hcp1200_processed/2021/Processed/${SUB}"
 
for (( c=$START; c<=$END; c++ ))
do
seedid="$(($c))"
seeddata=$(sed -n "${seedid}p" ${SEEDS})

tckgen ${PARENTDIR}/FOD.mif ${WORKDIR}/tracts_ALLGENEseeds/thal_seed_${c}.tck -seed_sphere ${seeddata},0 -select 5000 -seed_unidirectional -act ${PARENTDIR}/ACT.nii -backtrack

	for PARC in random500; do

		PARCPATH="${PARCDIR}/${PARC}_acpc.nii"
		tck2connectome ${WORKDIR}/tracts_ALLGENEseeds/thal_seed_${c}.tck $PARCPATH ${WORKDIR}/tracts_ALLGENEseeds/${PARC}/thal_seed_${c}_cp -vector -assignment_radial_search 5 -force

	done

done

fi
