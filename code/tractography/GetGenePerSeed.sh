#!/bin/env bash

module load mrtrix

DIR="/projects/kg98/stuarto/SeedReg"
MM=1.75

GENE_DIR="/projects/hf49/Stuart/AHBA_wholebrain"
GENE_LIST="${DIR}/GenesKept.txt"

ngene=$(wc -l ${GENE_LIST} | awk '{ print $1 }')

SEED_GENEDATA_DIR="${DIR}/seeds_${MM}mm_tracts_gene"

SEEDS="${DIR}/seeds_${MM}mm_tracts.tck"
echo "${ngene} genes"
mkdir ${SEED_GENEDATA_DIR}

for IND in $(seq 1 $ngene); do

GENE=$(sed -n "${IND}p" ${GENE_LIST})

tcksample -nointerp ${SEEDS} ${GENE_DIR}/${GENE}/${GENE}_mRNA.nii ${SEED_GENEDATA_DIR}/seeds_${MM}mm_tracts_${GENE}_gene.txt

done
