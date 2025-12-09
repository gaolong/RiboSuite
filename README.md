# RiboSuite

**RiboSuite** is a Nextflow DSL2–based, end-to-end pipeline for comprehensive **Ribo-seq data analysis**, including preprocessing, alignment, P-site assignment, 3-nt periodicity QC, ORF-level quantification, and downstream integration with complementary datasets (RNA-seq, proteomics, microprotein predictions).

RiboSuite is designed with **modularity**, **scalability**, and **high-performance computing (HPC)** environments in mind, following nf-core style conventions.

---

## Features

- End-to-end Ribo-seq processing using Nextflow DSL2
- Modular architecture compatible with nf-core patterns
- Supports alignment (STAR / Bowtie2)
- Automated P-site offset estimation
- 3-nt periodicity QC + metagene analysis
- ORF-level quantification and summarization
- Containerized execution (Conda, Mamba, Singularity)
- Scales from laptop to large HPC clusters
- Fully reproducible and version-controlled workflow

---

## Installation

You need:

- Nextflow (≥ 23.x)
- Mamba/Conda or Singularity
- A clone of this repository

Clone the repository:

```bash
git clone git@github.com:gaolong/RiboSuite.git
cd RiboSuite
