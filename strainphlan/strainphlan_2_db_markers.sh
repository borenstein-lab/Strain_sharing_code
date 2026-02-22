#!/bin/bash

# ==============================================================
# Create tree for all clades with strainphlan 
# ==============================================================

# Last update: 12/11/2025
# -- Usage example:
# bash strainphlan_2_db_markers.sh "t__SGB4303" 1 5               "t__SGB1557" "t__SGB14205" "t__SGB4303"

# -- Mandatory inputs
# CLADE - Sample name

# -- Arguments
CLADE="$1"

TOP_DIR="../../../strainphlan"	# Top directory for strainphlan analysis
MPA_DB_PATH="../../../../../PROGRAMS/conda/envs/mpa4/lib/python3.7/site-packages/metaphlan/metaphlan_databases/mpa_vJan21_CHOCOPhlAnSGB_202103.pkl"	
OVERRIDE=false	# Should existing files be replaced? If true - existing files will be overriden at any processing stage. If false - skip steps where output files already exist. 

# ==============================================================
# Setup 
# ==============================================================
if [ -z ${CLADE} ]; then echo "Missing argument: CLADE"; exit 1; fi
echo "CLADE: ${CLADE}"

# Activate Metaphlan4 environment
conda activate mpa4

# ==============================================================
# Extract DB markers for clade
# ==============================================================
DB_MARKERS_DIR="${TOP_DIR}/db_markers"

echo $(date) "Extract DB markers for clade"
mkdir -p ${DB_MARKERS_DIR}
cd ${DB_MARKERS_DIR}
if [[ ! -e "${DB_MARKERS_DIR}/${CLADE}.fna" ]]; then
    extract_markers.py -c ${CLADE} -o ${DB_MARKERS_DIR}/ -d ${MPA_DB_PATH}
else
    echo "DB markers for clade ${CLADE} already exist. Skipping extraction."
fi
# ==============================================================

conda deactivate
exit 0