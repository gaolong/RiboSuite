process FRAME_PERIODICITY_QC {

    tag "$sample_id"

    conda "bioconda::pysam=0.22.0 conda-forge::matplotlib conda-forge::pandas"

    publishDir "${params.outdir}/qc/frame_periodicity",
        mode: 'copy',
        saveAs: { file -> "${sample_id}/${file}" }

    input:
        tuple val(sample_id),
              path(bam),
              path(bai),
              path(psite_offset)
        path gtf

    output:
        tuple val(sample_id),
              path("${sample_id}.periodicity.by_region_by_length.tsv"),
              emit: periodicity_tsv

        path("${sample_id}.periodicity.by_region_by_length.heatmap.png"),
             optional: true,
             emit: periodicity_png

    script:
    """
    set -euo pipefail

    python ${moduleDir}/frame_periodicity.py \
        --bam ${bam} \
        --gtf ${gtf} \
        --offset ${psite_offset} \
        --sample ${sample_id}
    """
}
