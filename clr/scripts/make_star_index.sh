#!/bin/bash
#PBS -N star_genome_generate
#PBS -q super
#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=8:mem=40gb:ngpus=0
#PBS -k oe

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
BASE=/data/holwm/home/gaolong/projects/RiboSuite/reference/human/GRCh38

GENOME_DIR=${BASE}/star
FASTA=${BASE}/GRCh38.primary_assembly.genome.fa
GTF=${BASE}/gencode.v45.annotation.gtf

# -----------------------------
# Create STAR index directory
# -----------------------------
mkdir -p ${GENOME_DIR}

# -----------------------------
# Build STAR genome index
# -----------------------------
STAR --runThreadN 8 \
     --runMode genomeGenerate \
     --genomeDir ${GENOME_DIR} \
     --genomeFastaFiles ${FASTA} \
     --sjdbGTFfile ${GTF} \
     --sjdbOverhang 100

# -----------------------------
# Done
# -----------------------------
echo "STAR genomeGenerate finished successfully"

