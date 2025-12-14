nextflow.enable.dsl=2

process BOWTIE2_FILTER {
    tag "$sample_id"
    conda "bioconda::bowtie2=2.5.2 bioconda::samtools=1.18"

    input:
    tuple val(sample_id), path(reads)
    path index_prefix

    output:
    tuple val(sample_id), path("${sample_id}.clean.fastq.gz"), emit: clean_reads
    tuple val(sample_id), path("${sample_id}.contam.bam"), emit: contam_bam

    script:
    """
    # Align to contaminant index
    bowtie2 \
        -x $index_prefix \
        -U $reads \
        --very-sensitive \
        -p ${task.cpus} \
        2> ${sample_id}.bowtie2.log \
    | samtools view -b -F 4 -o ${sample_id}.contam.bam

    # Extract unmapped reads (clean reads)
    samtools fastq -f 4 ${sample_id}.contam.bam \
      | gzip > ${sample_id}.clean.fastq.gz
    """
}
