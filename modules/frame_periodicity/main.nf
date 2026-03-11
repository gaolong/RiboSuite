nextflow.enable.dsl = 2

process FRAME_PERIODICITY_QC {

    tag "${meta.sample_id}"

    conda "bioconda::pysam=0.22.0 conda-forge::matplotlib conda-forge::pandas"

    publishDir "${params.outdir}/ribo/qc/frame_periodicity",
        mode: 'copy',
        saveAs: { file -> "${meta.sample_id}/${file}" }

    input:
        tuple val(meta),
              path(bam),
              path(bai),
              path(psite_offset),
              path(gtf)

    output:
        tuple val(meta),
              path("${meta.sample_id}.periodicity.by_region_by_length.tsv"),
              emit: periodicity_tsv

        tuple val(meta),
              path("${meta.sample_id}.periodicity.by_region_by_length.heatmap.png"),
              optional: true,
              emit: periodicity_png

    script:
    """
    set -euo pipefail

    python ${moduleDir}/frame_periodicity.py \
        --bam ${bam} \
        --gtf ${gtf} \
        --psite_offsets ${psite_offset} \
        --sample ${meta.sample_id}
    """
}