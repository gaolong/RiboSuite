process PSITE_OFFSET {

    tag "$sample_id"

    conda "bioconda::pysam=0.22.1 python=3.10"

    input:
    tuple val(sample_id), path(bam)
    path gtf

    output:
    tuple val(sample_id), path("${sample_id}.psite_offsets.tsv")

    script:
    """
    python ${projectDir}/modules/psite_offset/psite_offset.py \
      --bam $bam \
      --gtf $gtf \
      --sample_id $sample_id \
      --out ${sample_id}.psite_offsets.tsv
    """
}