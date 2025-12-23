#!/bin/bash
#PBS -N contam_filter_test
#PBS -q workq
#PBS -l walltime=02:00:00
#PBS -l select=1:ncpus=4:mem=8gb
#PBS -k oe
#PBS -o /data/holwm/home/gaolong/projects/RiboSuite/clr/std/contam_filter_test.o
#PBS -e /data/holwm/home/gaolong/projects/RiboSuite/clr/std/contam_filter_test.e

cd $PBS_O_WORKDIR

source /data/holwm/home/gaolong/miniforge3/etc/profile.d/conda.sh
conda activate ribo-dev

IN_FASTQ=/data/holwm/home/gaolong/projects/RiboSuite/results/fastp_test/SRR1173905.trimmed.fastq.gz
IDX=/data/holwm/home/gaolong/projects/RiboSuite/reference/human/GRCh38/contam/bowtie2_idx/human_contam
OUTDIR=/data/holwm/home/gaolong/projects/RiboSuite/results/contam_filter_test

mkdir -p ${OUTDIR}

bowtie2 \
  -x ${IDX} \
  -U ${IN_FASTQ} \
  --threads 4 \
  --un-gz ${OUTDIR}/SRR1173905.noncontam.fastq.gz \
  -S ${OUTDIR}/SRR1173905.contam.sam