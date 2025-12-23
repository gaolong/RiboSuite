#!/bin/bash
#PBS -N fastp_test
#PBS -q workq
#PBS -l walltime=02:00:00
#PBS -l select=1:ncpus=4:mem=8gb:ngpus=0
#PBS -k oe
#PBS -o /data/holwm/home/gaolong/projects/RiboSuite/clr/std/fastp_test.o
#PBS -e /data/holwm/home/gaolong/projects/RiboSuite/clr/std/fastp_test.e

# -----------------------------
# Go to submission directory
# -----------------------------
cd $PBS_O_WORKDIR

# -----------------------------
# Activate conda environment
# -----------------------------
source ~/.bashrc
conda activate ribo-dev

# -----------------------------
# Define paths
# -----------------------------
FASTQ=/data/holwm/home/gaolong/projects/RiboSuite/testdata/human_ribo/SRR1173905.small.fastq.gz

OUTDIR=/data/holwm/home/gaolong/projects/RiboSuite/results/fastp_test
STD=/data/holwm/home/gaolong/projects/RiboSuite/clr/std

mkdir -p ${OUTDIR}
mkdir -p ${STD}

# -----------------------------
# Run fastp (single-end)
# -----------------------------
fastp \
  -i ${FASTQ} \
  -o ${OUTDIR}/SRR1173905.trimmed.fastq.gz \
  --detect_adapter_for_pe \
  --thread 4 \
  -h ${OUTDIR}/SRR1173905.fastp.html \
  -j ${OUTDIR}/SRR1173905.fastp.json

# -----------------------------
# Done
# -----------------------------
echo "fastp test finished successfully"
