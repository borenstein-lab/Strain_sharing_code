#!/bin/bash

# ==============================================================
# Create tree for all clades with strainphlan 
# ==============================================================

# Last update: 12/11/2025
# -- Usage example:
# bash strainphlan_3_tree_ngd.sh "t__SGB4303" 1 5               

# -- Mandatory inputs
# CLADE - Sample name
# MARKER_IN_N_SAMPLES - Minimum number of samples a marker must be present in to be included in the tree
# SAMPLE_WITH_N_MARKERS - Minimum number of markers a sample must have to be included in

# -- Arguments
CLADE="$1"
MARKER_IN_N_SAMPLES="$2" #1
SAMPLE_WITH_N_MARKERS="$3" #8	     # Minimum number of markers a sample must have to be included in the tree

TOP_DIR="../../..//strainphlan"	# Top directory for strainphlan analysis
BASE_DIR="../../../strainphlan/results/sample_with_n_markers_${SAMPLE_WITH_N_MARKERS}_marker_in_n_samples_${MARKER_IN_N_SAMPLES}"	# Base directory for processing
TREE_DIR="${BASE_DIR}/tree_output"		# Output folder. Will contain all processed results and intermediate files.
LONG_DIR="${BASE_DIR}/ngd_long"	# Output folder for long ngd files
MAT_DIR="${BASE_DIR}/ngd_mat"	# Output folder for mat ngd files
OVERRIDE=false	# Should existing files be replaced? If true - existing files will be overriden at any processing stage. If false - skip steps where output files already exist. 
N_THREADS=16	# Number of threads - used by any program that allows parallelization

# ==============================================================
# Setup 
# ==============================================================
if [ -z ${CLADE} ]; then echo "Missing argument: CLADE"; exit 1; fi
if [ -z ${MARKER_IN_N_SAMPLES} ]; then echo "Missing argument: MARKER_IN_N_SAMPLES"; exit 1; fi
if [ -z ${SAMPLE_WITH_N_MARKERS} ]; then echo "Missing argument: SAMPLE_WITH_N_MARKERS"; exit 1; fi

echo "CLADE: ${CLADE}"
echo "SAMPLE_WITH_N_MARKERS: ${SAMPLE_WITH_N_MARKERS}"
echo "MARKER_IN_N_SAMPLES: ${MARKER_IN_N_SAMPLES}"

# Activate Metaphlan4 environment
conda activate mpa4

mkdir -p ${BASE_DIR}

# ==============================================================
# Generation of best tree file
# ==============================================================

mkdir -p ${BASE_DIR}
mkdir -p ${TREE_DIR}
mkdir -p ${TREE_DIR}/${CLADE}
cd ${BASE_DIR}

if [[ ! -f "${TREE_DIR}/${CLADE}/RAxML_bestTree.${CLADE}.StrainPhlAn4.tre" ]]; then 
    echo "Working on tree for clade ${CLADE}"
    strainphlan -s ../../consensus_markers/*.pkl \
                -m ../../db_markers/${CLADE}.fna \
                -o ${TREE_DIR}/${CLADE} \
                -n ${N_THREADS} \
                -c ${CLADE} \
                --mutation_rates \
                --marker_in_n_samples ${MARKER_IN_N_SAMPLES} \
                --sample_with_n_markers ${SAMPLE_WITH_N_MARKERS} \
                --phylophlan_mode accurate
else
    echo "Tree ${CLADE} was already created / something went wrong."
fi

if [[ -f "${TREE_DIR}/${CLADE}/RAxML_bestTree.${CLADE}.StrainPhlAn4.tre" ]]; then
    echo "Tree saved in ${TREE_DIR}/${CLADE}/RAxML_bestTree.${CLADE}.StrainPhlAn4.tre"
else
    echo "Files within clade folder:"
    ls "${TREE_DIR}/${CLADE}"
    rm -rf "${TREE_DIR}/${CLADE}"
    echo "No output files were created, ${CLADE} was removed. Exiting."
    exit 1
fi

# ==============================================================
# Generation of pairwise distance files
# ==============================================================

mkdir -p "${LONG_DIR}"
mkdir -p "${MAT_DIR}"

cd ${TREE_DIR}/${CLADE}
if [[ ! -f "${LONG_DIR}/${CLADE}_nGD.tsv" ]]; then
    python ../../../../../code/strainphlan/tree_pairwisedists.py -n RAxML_bestTree.${CLADE}.StrainPhlAn4.tre ${LONG_DIR}/${CLADE}_nGD.tsv
    echo "Generated ngd long file for clade ${CLADE} saved to ${LONG_DIR}/${CLADE}_nGD.tsv"
else
    echo "Ngd long file for clade ${CLADE} already exists. Skipping generation."
fi

if [[ ! -f "${MAT_DIR}/${CLADE}_nGD_mat.tsv" ]]; then
    python ../../../../../code/strainphlan/tree_pairwisedists.py -n -m RAxML_bestTree.${CLADE}.StrainPhlAn4.tre ${MAT_DIR}/${CLADE}_nGD_mat.tsv
    echo "Generated ngd mat file for clade ${CLADE} saved to ${MAT_DIR}/${CLADE}_nGD_mat.tsv"
else
    echo "Ngd mat file for clade ${CLADE} already exists. Skipping"
fi

conda deactivate
exit 0