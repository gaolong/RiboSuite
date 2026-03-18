nextflow.enable.dsl = 2

process PSITE_SHIFT_BAM {

    tag "${meta.sample_id}"

    conda "conda-forge::python=3.11 conda-forge::pandas bioconda::pysam bioconda::samtools=1.18"

    publishDir "${params.outdir}/ribo/align/psite_shift_bam",
        mode: 'copy',
        saveAs: { file -> "${meta.sample_id}/${file}" }

    input:
        tuple val(meta), path(bam), path(bai), path(offset_tsv)

    output:
        tuple val(meta),
              path("${meta.sample_id}.psite_shifted.bam"),
              path("${meta.sample_id}.psite_shifted.bam.bai"),
              path("${meta.sample_id}.psite_shifted.offsets.tsv"),
              emit: shifted

    script:
    """
    set -euo pipefail

    python ${moduleDir}/psite_shift_bam.py \\
        --bam ${bam} \\
        --offsets ${offset_tsv} \\
        --sample ${meta.sample_id} \\
        --out_bam ${meta.sample_id}.psite_shifted.unsorted.bam \\
        --out_offsets ${meta.sample_id}.psite_shifted.offsets.tsv

    samtools sort -o ${meta.sample_id}.psite_shifted.bam ${meta.sample_id}.psite_shifted.unsorted.bam
    samtools index ${meta.sample_id}.psite_shifted.bam

    rm -f ${meta.sample_id}.psite_shifted.unsorted.bam
    """
}