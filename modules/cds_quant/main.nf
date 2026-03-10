nextflow.enable.dsl = 2

process CDS_QUANT {

    tag "${meta.sample_id}"

    conda "conda-forge::python=3.11 conda-forge::pandas bioconda::pysam"

    publishDir "${params.outdir}/ribo/cds_quant",
        mode: 'copy',
        pattern: "${meta.sample_id}.cds_quant*.tsv"

    input:
        tuple val(meta), path(bam), path(bai), path(gtf), path(offsets)

    output:
        tuple val(meta),
              path("${meta.sample_id}.cds_quant.gene.tsv"),
              path("${meta.sample_id}.cds_quant.stats.tsv")

    script:
    """
    set -euo pipefail

    python ${moduleDir}/cds_quant.py \
        --bam ${bam} \
        --gtf ${gtf} \
        --psite_offsets ${offsets} \
        --sample ${meta.sample_id} \
        --out_prefix ${meta.sample_id}.cds_quant
    """
}