#!/bin/bash
# Created by Michal Bukowski (michal.bukowski@tuta.io) under GPL-3.0 license

# Prepares complete test input data for rnaseq-pipeline-2 pipeline by fetching
# reads from 12 Illumina MiSeq runs from NCBI SRA to reads/ directory.
# These archives encompass 2 wild type (wt) vs. mutant (mt) comparisons,
# in logarithmic (lg) and late (lt) growth phases, each in 3 replicas.

set -euo pipefail

# SRA run accession numbers for the read archives and final file names, which
# follow the pattern: strain_genotype_growth-phase_replica_RN.fastq.gz,
# where N is 1 or 2 and refers to reads 1 and reads 2 from read pairs.
declare -A SRA_RUNACC
SRA_RUNACC=(
    ["SRR17659943"]="rn_wt_lg_1_R"
    ["SRR22473284"]="rn_wt_lg_2_R"
    ["SRR22473283"]="rn_wt_lg_3_R"
    ["SRR17659941"]="rn_mt_lg_1_R"
    ["SRR22473279"]="rn_mt_lg_2_R"
    ["SRR22473280"]="rn_mt_lg_3_R"
    ["SRR17659942"]="rn_wt_lt_1_R"
    ["SRR22473287"]="rn_wt_lt_2_R"
    ["SRR22473288"]="rn_wt_lt_3_R"
    ["SRR17659940"]="rn_mt_lt_1_R"
    ["SRR22473282"]="rn_mt_lt_2_R"
    ["SRR22473281"]="rn_mt_lt_3_R"
)

# Check the path for saving reads to. If fails, exit with code 1.
if [[ -e "reads" ]]; then
    echo "Cannot create reads/ output directory. The path exists." \
         "Rename of remove the directory and start again." >&2
    exit 1
fi

# Warn the user about the size of dowloaded data. On rejection, exit with code 0.
echo "This script will fetch 12 read archives to reads/ directory (~3.7 GB)."
echo "It may take a while..."
answer=""
while [[ "${answer}" != "yes" && "${answer}" != "no" ]]; do
    echo "Do you want to proceed? (yes/no)"
    read -r answer
done
if [[ "${answer}" == "no" ]]; then
    echo "Exiting..."
    exit 0
fi

# Download test read archives. If fails at any point, exit with code 1.
count="${#SRA_RUNACC[@]}"
echo "Fetching ${count} test read archives..."
mkdir -p "reads"
fastq-dump -v --gzip --split-files --outdir "reads" ${!SRA_RUNACC[@]}
for runacc in "${!SRA_RUNACC[@]}"; do
    fname="${SRA_RUNACC[$runacc]}"
    for n in 1 2; do
        mv "reads/${runacc}_${n}.fastq.gz" "reads/${fname}${n}.fastq.gz"
    done
done
echo "Done"

