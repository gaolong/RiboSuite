process CDS_QUANT {

    tag "$sample_id"

    conda "conda-forge::python=3.11 conda-forge::pandas bioconda::pysam"

    publishDir "${params.outdir}/ribo/cds_quant",
        mode: 'copy',
        pattern: "${sample_id}.cds_quant*.tsv"

    input:
        tuple val(sample_id), path(bam), path(bai)
        path gtf
        path offsets

    output:
        tuple val(sample_id),
              path("${sample_id}.cds_quant.gene.tsv"),
              path("${sample_id}.cds_quant.stats.tsv")

    script:
    """
    set -euo pipefail

    python ${moduleDir}/cds_quant.py \
        --bam $bam \
        --gtf $gtf \
        --offsets $offsets \
        --sample $sample_id \
        --out_prefix ${sample_id}.cds_quant
    """
}
