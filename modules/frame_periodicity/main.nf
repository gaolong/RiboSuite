process FRAME_PERIODICITY_QC {

    tag "$sample_id"

    conda "bioconda::pysam=0.22.0 conda-forge::matplotlib conda-forge::pandas"

    input:
        tuple val(sample_id),
              path(bam),
              path(bai),
              path(psite_offset)
        path gtf

    output:
        tuple val(sample_id),
              path("${sample_id}.frame_counts.tsv"),
              path("${sample_id}.frame_fraction.tsv"),
              path("${sample_id}.frame_periodicity.png")

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
