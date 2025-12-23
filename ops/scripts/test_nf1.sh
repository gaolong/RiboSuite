#!/bin/bash
#PBS -N ribosuite
#PBS -q workq
#PBS -l select=1:ncpus=8:mem=40gb
#PBS -l walltime=12:00:00
#PBS -j oe
#PBS -o /data/holwm/home/gaolong/projects/RiboSuite/ops/logs/ribosuite.o

cd /data/holwm/home/gaolong/projects/RiboSuite

source ~/.bashrc
conda activate ribo-dev

nextflow run . \
  --reads "testdata/human_ribo/SRR1173905.small.fastq.gz" \
  -profile pbs 
