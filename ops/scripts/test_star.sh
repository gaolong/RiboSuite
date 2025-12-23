#!/bin/bash
#PBS -N star_align_test
#PBS -q workq
#PBS -l walltime=02:00:00
#PBS -l select=1:ncpus=8:mem=32gb
#PBS -k oe
#PBS -o /data/holwm/home/gaolong/projects/RiboSuite/ops/std/star_align_test.o
#PBS -e /data/holwm/home/gaolong/projects/RiboSuite/ops/std/star_align_test.e

cd $PBS_O_WORKDIR

source /data/holwm/home/gaolong/miniforge3/etc/profile.d/conda.sh
conda activate ribo-dev

GENOME_DIR=/data/holwm/home/gaolong/projects/RiboSuite/reference/human/GRCh38/star
IN_FASTQ=/data/holwm/home/gaolong/projects/RiboSuite/results/contam_filter_test/SRR1173905.noncontam.fastq.gz
OUTDIR=/data/holwm/home/gaolong/projects/RiboSuite/results/star_test

mkdir -p ${OUTDIR}

STAR \
  --genomeDir ${GENOME_DIR} \
  --readFilesIn ${IN_FASTQ} \
  --readFilesCommand zcat \
  --runThreadN 8 \
  --outFileNamePrefix ${OUTDIR}/SRR1173905. \
  --outSAMtype BAM SortedByCoordinate \
  --outFilterMismatchNmax 2 \
  --alignEndsType EndToEnd

echo "STAR alignment finished"
