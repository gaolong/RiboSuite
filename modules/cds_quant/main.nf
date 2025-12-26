process CDS_QUANT {

    tag "$sample_id"

    conda "bioconda::pysam bioconda::pandas"

    input:
    tuple val(sample_id), path(bam), path(bai)
    path gtf
    path offsets

    output:
    tuple val(sample_id),
          path("${sample_id}.cds_quant.transcript.tsv"),
          path("${sample_id}.cds_quant.gene.tsv")

    script:
    """
    cds_quant.py \
        --bam $bam \
        --gtf $gtf \
        --offsets $offsets \
        --sample $sample_id \
        --out_prefix ${sample_id}.cds_quant
    """
}
