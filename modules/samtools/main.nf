nextflow.enable.dsl=2

process SAMTOOLS_SORT_INDEX {
    tag "$sample_id"
    conda "bioconda::samtools=1.18"

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id),
          path("${sample_id}.sorted.bam"),
          path("${sample_id}.sorted.bam.bai")

    script:
    """
    samtools sort -@ ${task.cpus} -o ${sample_id}.sorted.bam $bam
    samtools index ${sample_id}.sorted.bam
    """
}