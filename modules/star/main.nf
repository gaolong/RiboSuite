nextflow.enable.dsl=2

process STAR_ALIGN {
    tag "$sample_id"
    conda "bioconda::star=2.7.11b"

    input:
    tuple val(sample_id), path(reads)
    path star_index

    output:
    tuple val(sample_id), path("${sample_id}.Aligned.out.bam")

    script:
    """
    STAR \
      --genomeDir $star_index \
      --readFilesIn $reads \
      --readFilesCommand zcat \
      --outSAMtype BAM Unsorted \
      --outFileNamePrefix ${sample_id}. \
      --runThreadN ${task.cpus}
    """
}