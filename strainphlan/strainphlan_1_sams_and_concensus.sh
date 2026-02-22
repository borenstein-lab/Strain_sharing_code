#!/bin/bash

# ==============================================================
# Process a single metaphlan sample with median files 
# ==============================================================

# Last update: 26/11/2025

# The following folders are created by default inside the output folder:
# `logs` for all log files created during processing
# `sams` 
# `profiles`
# `bowtie`
# `consensus_markers`

# Usage example:
#  bash strainphlan_1_sams_and_concensus.sh "RM-125-vaginal_CP04600_L001" > log_1066422.txt

# ==============================================================
# Arguments
# ==============================================================

# Mandatory command line arguments
# -------------------------------------------
# SAMPLE_ID - Sample name, expected to have 2 associated files in the input folder, namely <SAMPLE_ID>_<FASTQ_FWD_SUFFIX> and <SAMPLE_ID>_<FASTQ_REV_SUFFIX>

# Other configurations with default settings
# -------------------------------------------
# Note: These can be optionally overriden by passing them as additional command line arguments (KEY=VALUE)

# -- General
IN_DIR=../../../strainphlan/ # fastqs_post_qc		# Input folder (contains fastq files after QC)
OUT_DIR=../../../strainphlan		# Output folder. Will contain all processed results and intermediate files.
SAMPLE_ID="$1"
OVERRIDE=false	# Should existing files be replaced? If true - existing files will be overriden at any processing stage. If false - skip steps where output files already exist. 
N_THREADS=16	# Number of threads - used by any program that allows parallelization
FASTQ_CLEAN_FINAL="${IN_DIR}/${SAMPLE_ID}_CLEAN.fastq.gz"

SAMS=true
CONCENSUS=true
# ==============================================================
# Setup 
# ==============================================================
# Activate Metaphlan4 environment
conda activate mpa4

if [[ "${SAMS}" = true ]]; then
    echo $(date) "Creating sams files for SAMPLE_ID: $SAMPLE_ID"
    mkdir -p "${OUT_DIR}/sams"
    mkdir -p "${OUT_DIR}/bowtie2"
    mkdir -p "${OUT_DIR}/profiles"
fi

# Start
echo $(date) "Processing SAMPLE_ID: $SAMPLE_ID"

# ==============================================================
# Metaphlan4
# ==============================================================

if [[ "${SAMS}" = true ]]; then
    echo $(date) "Running Metaphlan4 to create sams files"
    metaphlan $FASTQ_CLEAN_FINAL \
        --input_type fastq \
        -o ${OUT_DIR}/profiles/${SAMPLE_ID}_metaphlan4_bugs_list.tsv \
        -s ${OUT_DIR}/sams/${SAMPLE_ID}.sam.bz2 \
        --bowtie2out ${OUT_DIR}/bowtie2/${SAMPLE_ID}.bowtie2.bz2 \
        --nproc $N_THREADS \
        -t rel_ab_w_read_stats \
        --unclassified_estimation \
        --no_map

    if [[ -e "${OUT_DIR}/sams/${SAMPLE_ID}.sam.bz2" ]]; then
      echo $(date) "Metaphlan4 completed successfully."
    else
      echo $(date) "Metaphlan4 failed to produce the expected output file."
    fi
fi

# ==============================================================
# Concensus markers extraction
# ==============================================================

if [[ "${CONCENSUS}" = true ]]; then
  mkdir -p ${OUT_DIR}/consensus_markers
  echo $(date) "Extracting concensus markers"
  consensus_file_path="${OUT_DIR}/consensus_markers/${SAMPLE_ID}.pkl"
  if [[-e "${consensus_file_path}" ]]; then
    echo "File exists."
  else
    sample2markers.py -i "${IN_DIR}/sams/${SAMPLE_ID}.sam.bz2" -o "${OUT_DIR}/consensus_markers" -n "${N_THREADS}"
    echo $(date) "Completed concensus markers extraction for SAMPLE_ID: $SAMPLE_ID"
  fi
fi

# ==============================================================
conda deactivate
echo $(date) "Completed StrainPhlAn first steps for SAMPLE_ID: $SAMPLE_ID"
exit 0