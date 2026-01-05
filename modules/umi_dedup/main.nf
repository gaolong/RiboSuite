process UMI_DEDUP {

    tag "$sample_id"

    conda "bioconda::umi_tools=1.1.4 bioconda::samtools=1.20"

    publishDir "${params.outdir}/align",
        mode: 'copy',
        pattern: "${sample_id}.umi_dedup.{bam,bam.bai,log}"
        

    input:
    tuple val(sample_id), path(bam), path(bai)

    output:
    tuple val(sample_id),
          path("${sample_id}.umi_dedup.bam"),
          path("${sample_id}.umi_dedup.bam.bai"),
          path("${sample_id}.umi_dedup.log")

    script:
    """
    set -euo pipefail

    umi_tools dedup \
        --stdin ${bam} \
        --stdout ${sample_id}.umi_dedup.bam \
        --log ${sample_id}.umi_dedup.log

    samtools index ${sample_id}.umi_dedup.bam
    """
}
