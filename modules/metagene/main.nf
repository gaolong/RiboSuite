process METAGENE_QC {

    tag "$sample_id"
    label 'small'

    conda "bioconda::pysam=0.22.0 bioconda::htslib=1.20 bioconda::samtools=1.20 \
           conda-forge::python=3.10 conda-forge::pandas=2.2.2 \
           conda-forge::numpy=2.0.1 conda-forge::matplotlib=3.9.1"

    publishDir "${params.outdir}/metagene",
        mode: 'copy',
        pattern: "*.{png,tsv}"

    input:
        tuple val(sample_id),
              path(bam),
              path(bai),
              path(offset)
        path gtf

    output:
        tuple val(sample_id),
              path("${sample_id}.metagene.start.tsv"),
              path("${sample_id}.metagene.stop.tsv"),
              path("${sample_id}.metagene.png")

    script:
    """
    set -euo pipefail

    python ${moduleDir}/metagene.py \
      --bam ${bam} \
      --gtf ${gtf} \
      --offset ${offset} \
      --sample ${sample_id} \
      --window ${params.metagene_window ?: 50} \
      --max_transcripts_per_gene ${params.metagene_max_tx_per_gene ?: 1}
    """
}
